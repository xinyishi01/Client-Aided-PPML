#include "utils.hpp"
#include "Layer.hpp"


class Layer_Server : public Layer {
public:

	Mat x_mask, W_mask, Go_mask;
	Mat forward_tp, Wg_tp, backward_tp;
	
	Layer_Server() {}
	
	Layer_Server(int _party, int _D_in, int _D_out, string _id, long int _lr, PRF *_prf_server) {
		party = _party;
		D_in = _D_in;
		D_out = _D_out;
		id = _id;
		lr = _lr;
		W.resize(D_in, D_out);
		W_grad.resize(D_in, D_out);
		bias.resize(1, D_out);
		bias_grad.resize(1, D_out);
		z1.resize(B, D_out);
		z2.resize(B, D_out);
		b1.resize(B, D_out);
		tau1.resize(B, D_out);
		tau2.resize(B, D_out);
		select_comparison.resize(B, D_out);
		for (int i = 0; i < D_in; ++i) for (int j = 0; j < D_out; ++j) W(i, j) = myrand();
		for (int i = 0; i < D_out; ++i) bias(0, i) = 0;//myrand2();
		it = 0;
		
		x_mask.resize(B, D_in);
		W_mask.resize(D_in, D_out);
		Go_mask.resize(B, D_out);
		forward_tp.resize(B, D_out);
		Wg_tp.resize(D_in, D_out);
		backward_tp.resize(B, D_in);

		for (int i = L; i >= 0; i--)
			b_mask[i].resize(B, D_out);
		
		prf_server = _prf_server;
		
		act_time = mul_time = off_time = 0.;
		
		COMM = COMM_client = 0.;
	}


	void pre(PRF *_prf_client, PRF *_prf_client_same) {
		int t1 = mytime.start();
		prf_client = _prf_client;
		prf_client_same = _prf_client_same;
	
		cal_mat2(prf_client, x_mask);
		prf_client_same->cal(W_mask);
		cal_mat2(prf_client, Go_mask);
		
		cal_mat2(prf_client, z1);
		if(id != "L3") cal_mat2(prf_client, z2);

		if (id == "L3") {
			for (int i = L; i >= 0; i--)
				cal_mat2(prf_client, b_mask[i]);
		}
		
		COMM_client += merge_mat(forward_tp, prf_client, io_client, party);
		COMM_client += merge_mat(backward_tp, prf_client, io_client, party);
		
		Mat *temp = new Mat(B, D_in * D_out);
		
		COMM_client += merge_mat(*temp, prf_client, io_client, party);
		
		for (int i = 0, cnt = 0; i < D_in; ++i) for (int j = 0; j < D_out; ++j, ++cnt) {
			Wg_tp(i, j) = 0;
			for (int k = 0; k < B; ++k) Wg_tp(i, j) += (*temp)(k, cnt);
		}
		
		delete temp;
	
		exchange_comparisons(select_comparison);
	
		off_time += mytime.end(t1);
	}

	Mat forward(Mat &x) {
		in_x = x;
		int t1 = mytime.start();
		x_mask = reconstruct(in_x - x_mask, party, io_server); 
		COMM += 2 * in_x.rows() * in_x.cols() * sizeof(unsigned long int);
		W_mask = reconstruct(W - W_mask, party, io_server);
		COMM += 2 * W.rows() * W.cols() * sizeof(unsigned long int);
		Mat inner = inner_product(in_x, W, x_mask, W_mask, forward_tp);
		truncation_mat(inner);
		for (int i = 0; i < B; ++i) inner.row(i) = inner.row(i) + bias;
		mul_time += mytime.end(t1);
		
		t1 = mytime.start();
		Mat rs = activation_function_server(inner);
		act_time += mytime.end(t1);
		
		if(DEBUG_FORWARD) {
			Mat tmpx = reconstruct(in_x, party, io_server);
			Mat tmpW = reconstruct(W, party, io_server);
			Mat tmpb = reconstruct(bias, party, io_server);
			Mat predict_rs = plain_ReLU(tmpx, tmpW, tmpb);
			Mat rec_rs = reconstruct(rs, party, io_server);
			for(int i = 0; i < rec_rs.rows(); ++i) for (int j = 0; j < rec_rs.cols(); ++j) {
				if(abs((long int)rec_rs(i, j) - (long int)predict_rs(i, j)) > 2) {
					cout << id << "\n";
					cout << "(" << i << "," << j << "): \n" 
						<< (long int) predict_rs(i, j) << " " << (double) ((long int) predict_rs(i, j)) / (double)pow(2, L) << "\n"
						<< (long int)rec_rs(i, j) << " " << (double) ((long int) rec_rs(i, j)) / (double)pow(2, L) << "\n";
					
					assert(0);
				}
			}
			cerr << "ok\n";
		}
		
		return rs;
	}
	
	Mat backward(Mat& x) {
		++it;
		int t1 = mytime.start();
		Mat G_out(B, D_out);
		if(id != "L3") {
			G_out = cwise_product(x, z2, b1, party) + tau2;
			for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) if(select_comparison(i, j)) 
				G_out(i, j) = x(i, j) - G_out(i, j);
		}
		else G_out = x;
		act_time += mytime.end(t1);
		if(DEBUG_BACKWARD) {
			Mat rec_G_out = reconstruct(G_out, party, io_server);
			if(party == ALICE) debug_mat(rec_G_out, id + "_G_out.txt");
		}
		t1 = mytime.start();
		Go_mask = reconstruct(G_out - Go_mask, party, io_server);
		COMM += 2 * G_out.rows() * G_out.cols() * sizeof(unsigned long int);
		W_grad = inner_product(in_x.transpose(), G_out, x_mask.transpose(), Go_mask, Wg_tp);
		truncation_mat(W_grad);
		bias_grad = G_out.colwise().sum();
		Mat rs = inner_product(G_out, W.transpose(), Go_mask, W_mask.transpose(), backward_tp);
		truncation_mat(rs);
		mul_time += mytime.end(t1);
		return rs;
	}
	
	Mat softmax(Mat& X, Mat& Y) {
		//cout << "soft max...\n";
        int rows = B;
		int columns = K;
		for (int i = 0; i < rows; ++i) for (int j = 0; j < columns; ++j) X(i, j) = X(i, j) + 1;
		Mat tmp_ = X;
		Mat tmp__= X.rowwise().sum();
		Mat s(rows, columns), las(rows, columns), ans(rows, columns);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++) {
				s(i, j) = tmp__(i, 0);
				las(i, j) = tmp_(i, j);
				ans(i, j) = 0;
			}

		Mat z(rows, columns), v(rows, columns), b(rows, columns);
		for (int i = L; i >= 0; i--) {
			z = las - s;
			Mat v = sign_check(z, b, s, b_mask[i]);

			ans += (1UL << i) * b;
			las = 2 * (las - v);
		}
		
		if(false) {
			Mat rec_X = reconstruct(X, party, io_server);
			Mat rec_X_sum = reconstruct(tmp__, party, io_server);
			Mat rec_ans = reconstruct(ans, party, io_server);
			Mat cor_ans(B, K);
			for (int i = 0; i < B; ++i) for (int j = 0; j < K; ++j) {
				if(rec_X_sum(i, 0)) {
					cor_ans(i, j) = (rec_X(i, j) << L) / rec_X_sum(i, 0);
					if(cor_ans(i, j) != rec_ans(i, j)) {
						cout << i << "," << j << ": " << cor_ans(i, j) << " " << rec_ans(i, j) << "\n";
						assert(0);
					}
				}
			}
		}
		//cout << "done\n";
		return ans - Y;
	}

};
