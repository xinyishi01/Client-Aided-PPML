#include "utils.hpp"
#include "Layer.hpp"


class Layer1_Server : public Layer {
public:
	
	Mat z, tau, backward_z, backward_tau;
	
	Layer1_Server() {}
	
	Layer1_Server(int _party, int _D_in, int _D_out, string _id, long int _lr, PRF *_prf_server) {
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
		assert(D_in == D);
		it = 0;
		z.resize(D_in, D_out);
		tau.resize(B, D_out);
		backward_z.resize(B, D_out);
		backward_tau.resize(D_in, D_out);
		
		prf_server = _prf_server;
		
		act_time = mul_time = off_time = 0.;
		
		COMM = COMM_client = 0.;
	}

	void pre(PRF *_prf_client, PRF *_prf_client_same, Mat _in_x) {
		int t1 = mytime.start();
		prf_client = _prf_client;
		prf_client_same = _prf_client_same;
		in_x = _in_x;
	
		prf_client_same->cal(z);
		cal_mat2(prf_client, z1);
		cal_mat2(prf_client, z2);
		cal_mat2(prf_client, backward_z);
		
		COMM_client += merge_mat(tau, prf_client, io_client, party);
		
		Mat *temp = new Mat(B, D_in * D_out), *temp2 = new Mat(1, D_in * D_out);
		COMM_client += merge_mat(*temp, prf_client, io_client, party);
		(*temp2) = (*temp).colwise().sum();
		for (int i = 0, cnt = 0; i < D_in; ++i) for (int j = 0; j < D_out; ++j, ++cnt) {
			backward_tau(i, j) = (*temp2)(0, cnt);
		}
		delete temp;
		delete temp2;
		
		exchange_comparisons(select_comparison);
		off_time += mytime.end(t1);
	}
	Mat forward() {
		int t1 = mytime.start();
		Mat inner = inner_product_input(z, in_x, tau);
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
	
	void backward(Mat& x) {
		++it;
		Mat G_out(B, D_out);
		int t1 = mytime.start();
		G_out = cwise_product(x, z2, b1, party) + tau2;
		for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) if(select_comparison(i, j)) 
			G_out(i, j) = x(i, j) - G_out(i, j);
		act_time += mytime.end(t1);
		if(DEBUG_BACKWARD) {
			Mat rec_G_out = reconstruct(G_out, party, io_server);
			if(party == ALICE) debug_mat(rec_G_out, id + "_G_out.txt");
		}
		t1 = mytime.start();
		Mat alpha = reconstruct(G_out - backward_z, party, io_server);
		COMM += 2 * G_out.rows() * G_out.cols() * sizeof(unsigned long int);
		
		W_grad = (in_x.transpose() * alpha) + backward_tau;
		truncation_mat(W_grad);
		
		bias_grad = G_out.colwise().sum();
		mul_time += mytime.end(t1);
	}

};
