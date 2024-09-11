#include "utils.hpp"
#include "Layer.hpp"


class Layer_Server : public Layer {
public:

	Mat x_mask, W_mask, Go_mask;
	Mat forward_tp, Wg_tp, backward_tp;

/*--------------------------Malicious--------------------------*/
	Mat oX[rN], oW[rN], oTP[rN], oG[rN], obTP[rN], oWg[rN];
/*--------------------------Malicious--------------------------*/

	
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
		COMM = COMM_client = 0.;
		
		/*--------------------------Malicious--------------------------*/
		oNum = 0;
		dNum = 0;
		/*--------------------------Malicious--------------------------*/
	}

/*--------------------------Malicious--------------------------*/

bool cherk_malicious() {
    PRF *prf = new PRF;
    COMM_client += send_keys(prf, io_client, 1);
    
    Mat *nwX = new Mat[rM];
    Mat *nwW = new Mat[rM];
    Mat *nwTP = new Mat[rM];
    Mat *nwG = new Mat[rM];
    Mat *nwbTP = new Mat[rM];
    Mat *nwWg = new Mat[rM];
    
    for (int i = 0; i < rM; ++i) {
        nwX[i].resize(B, D_in); 
        nwW[i].resize(D_in, D_out);
        nwTP[i].resize(B, D_out);
        nwG[i].resize(B, D_out);
        nwbTP[i].resize(B, D_in);
        nwWg[i].resize(D_in, D_out);
        prf->cal(nwX[i]);
        prf->cal(nwW[i]);
        prf->cal(nwG[i]);
    }

	Mat temp0(rM, B * D_out), temp1(rM, B * D_in), temp2(rM, D_in * D_out);
	COMM_client += merge_mat_same(temp0, prf, io_client, party);
	COMM_client += merge_mat_same(temp1, prf, io_client, party);
	COMM_client += merge_mat_same(temp2, prf, io_client, party);
	
	
	for (int i = 0; i < rM; ++i) {
		unfold_mat(nwTP[i], temp0.row(i));
		unfold_mat(nwbTP[i], temp1.row(i));
		unfold_mat(nwWg[i], temp2.row(i));
	}
	
	delete prf;
    
    int id[rM];
    my_shuffle(id, rM);
    bool fl = 1;
    
    uint64 ans = 0, ans_;

	Mat openX(nwX[0].rows(), nwX[0].cols());
	Mat openW(nwW[0].rows(), nwW[0].cols());
	Mat openG(nwG[0].rows(), nwG[0].cols());
    for (int i = rk * rN; i < rM && fl; i++) {
        int o = id[i]; /*
        Mat openX = reconstruct(nwX[o], party, io_server);
        Mat openW = reconstruct(nwW[o], party, io_server);
        Mat openTP = reconstruct(nwTP[o], party, io_server);
        Mat openG = reconstruct(nwG[o], party, io_server);
        Mat openbTP = reconstruct(nwbTP[o], party, io_server);
        Mat openWg = reconstruct(nwWg[o], party, io_server);
        if((openX * openW) != openTP) fl = 0;
        if((openX.transpose() * openG) != openWg) fl = 0;
        if((openG * openW.transpose()) != openbTP) fl = 0;*/
        if(party == ALICE) {
        	COMM += receive_mat(openX, io_server);
        	COMM += receive_mat(openW, io_server);
        	COMM += receive_mat(openG, io_server);
        	openX = openX + nwX[o];
        	openW = openW + nwW[o];
        	openG = openG + nwG[o];
        	ans += linear_com(openX * openW - nwTP[o]);
        	ans += linear_com(openX.transpose() * openG - nwWg[o]);
        	ans += linear_com(openG * openW.transpose() - nwbTP[o]);
        }
        else {
        	COMM += send_mat(nwX[o], io_server);
        	COMM += send_mat(nwW[o], io_server);
        	COMM += send_mat(nwG[o], io_server);
        	ans += linear_com(-nwTP[o]);
        	ans += linear_com(-nwWg[o]);
        	ans += linear_com(-nwbTP[o]);
        }
    }
	
	//if(fl) {
	int tempsz = rN * (rk - 1);
	Mat *tempX = new Mat[tempsz];
	Mat *tempW = new Mat[tempsz];
	Mat *tempG = new Mat[tempsz];
	Mat *tempX_ = new Mat[tempsz];
	Mat *tempW_ = new Mat[tempsz];
	Mat *tempG_ = new Mat[tempsz];
	for (int i = 0, cnt = 0; i < rN; i++) {
		int e = id[i * rk];
		for (int j = 1; j < rk; j++, ++cnt) {
			int o = id[i * rk + j];
			tempX[cnt] = nwX[e] - nwX[o];
			tempW[cnt] = nwW[e] - nwW[o];
			tempG[cnt] = nwG[e] - nwG[o];
			tempX_[cnt].resize(tempX[cnt].rows(), tempX[cnt].cols());
			tempW_[cnt].resize(tempW[cnt].rows(), tempW[cnt].cols());
			tempG_[cnt].resize(tempG[cnt].rows(), tempG[cnt].cols());
		}
	}
	
	if(party == ALICE) {
		for (int i = 0; i < tempsz; ++i) {
			COMM += send_mat(tempX[i], io_server);
			COMM += send_mat(tempW[i], io_server);
			COMM += send_mat(tempG[i], io_server);
		}
		for (int i = 0; i < tempsz; ++i) {
			COMM += receive_mat(tempX_[i], io_server);
			COMM += receive_mat(tempW_[i], io_server);
			COMM += receive_mat(tempG_[i], io_server);
		}
	}
	else {
		for (int i = 0; i < tempsz; ++i) {
			COMM += receive_mat(tempX_[i], io_server);
			COMM += receive_mat(tempW_[i], io_server);
			COMM += receive_mat(tempG_[i], io_server);
		}
		for (int i = 0; i < tempsz; ++i) {
			COMM += send_mat(tempX[i], io_server);
			COMM += send_mat(tempW[i], io_server);
			COMM += send_mat(tempG[i], io_server);
		}
	}
	
	for (int i = 0; i < tempsz; ++i) {
		tempX[i] = tempX[i] + tempX_[i];
		tempW[i] = tempW[i] + tempW_[i];
		tempG[i] = tempG[i] + tempG_[i];
	}
	
		Mat rs(rN * (rk-1), 3);
		for (int i = 0, cnt = 0; i < rN; i++) {
			int e = id[i * rk];
		    for (int j = 1; j < rk; j++, ++cnt) {
		    	int o = id[i * rk + j];/*
		    	Mat tempX = reconstruct(nwX[e] - nwX[o], party, io_server);
		    	Mat tempW = reconstruct(nwW[e] - nwW[o], party, io_server);
		    	Mat tempG = reconstruct(nwG[e] - nwG[o], party, io_server);*/
		        rs(i * (rk - 1) + j - 1, 0) = linear_com(inner_product(nwX[e], nwW[e], tempX[cnt], tempW[cnt], nwTP[o]) - nwTP[e]);
		        rs(i * (rk - 1) + j - 1, 1) = linear_com(inner_product(nwX[e].transpose(), nwG[e], tempX[cnt].transpose(), tempG[cnt], nwWg[o]) - nwWg[e]);
		        rs(i * (rk - 1) + j - 1, 2) = linear_com(inner_product(nwG[e], nwW[e].transpose(), tempG[cnt], tempW[cnt].transpose(), nwbTP[o]) - nwbTP[e]);
		    }
		    oX[i] = nwX[e];
		    oW[i] = nwW[e];
		    oG[i] = nwG[e];
		    oTP[i] = nwTP[e];
		    oWg[i] = nwWg[e];
		    obTP[i] = nwbTP[e];
		}
		ans += linear_com(rs);
		
	delete[] tempX;
	delete[] tempW;
	delete[] tempG;
	delete[] tempX_;
	delete[] tempW_;
	delete[] tempG_;

		if(party == ALICE) {
			io_server->recv_data(&ans_, sizeof(uint64));
			io_server->send_data(&ans, sizeof(uint64));
			io_server->flush();
		}
		else {
			io_server->send_data(&ans, sizeof(uint64));
			io_server->flush();
			io_server->recv_data(&ans_, sizeof(uint64));
		}
		if(ans != -ans_) fl = 0;
	//}
	
    delete[] nwX;
    delete[] nwW;
    delete[] nwTP;
    delete[] nwG;
    delete[] nwbTP;
    delete[] nwWg;
    
    
    return fl;
}
		
/*--------------------------Malicious--------------------------*/	


	void pre(PRF *_prf_client, PRF *_prf_client_same) {
	/*--------------------------Malicious--------------------------*/
		
		
		if(oNum == 0) {
			//cout << id << " is generating triples...\n";
			if(!cherk_malicious()) {
				cout << "have found malicious!\n";
				assert(0);
			}
			oNum = 0;
			//cout << "OK.\n";
		}
		
		x_mask = oX[oNum]; //B, D_in
		W_mask = oW[oNum]; //D_in, D_out
		forward_tp = oTP[oNum]; //B, D_out
		Go_mask = oG[oNum];//B, D_out
		backward_tp = obTP[oNum];//B, D_in
		Wg_tp = oWg[oNum];//D_in, D_out
		
		if((++oNum) == rN) oNum = 0;
		/*--------------------------Malicious--------------------------*/
		prf_client = _prf_client;
		prf_client_same = _prf_client_same;
		
		cal_mat2(prf_client, z1);
		if(id != "L3") cal_mat2(prf_client, z2);
		
		if (id == "L3") {
			for (int i = L; i >= 0; i--)
				cal_mat2(prf_client, b_mask[i]);
		}
	
		exchange_comparisons(select_comparison);
	
	}

	Mat forward(Mat &x) {
		in_x = x;
		x_mask = reconstruct(in_x - x_mask, party, io_server); 
		COMM += 2 * in_x.rows() * in_x.cols() * sizeof(unsigned long int);
		W_mask = reconstruct(W - W_mask, party, io_server);
		COMM += 2 * W.rows() * W.cols() * sizeof(unsigned long int);
		Mat inner = inner_product(in_x, W, x_mask, W_mask, forward_tp);
		truncation_mat(inner);
		for (int i = 0; i < B; ++i) inner.row(i) = inner.row(i) + bias;
		Mat rs = activation_function_server(inner);
		
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
		
		Mat G_out(B, D_out);
		if(id != "L3") {
			G_out = cwise_product(x, z2, b1, party) + tau2;
			for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) if(select_comparison(i, j)) 
				G_out(i, j) = x(i, j) - G_out(i, j);
		}
		else G_out = x;
		if(DEBUG_BACKWARD) {
			Mat rec_G_out = reconstruct(G_out, party, io_server);
			if(party == ALICE) debug_mat(rec_G_out, id + "_G_out.txt");
		}
		Go_mask = reconstruct(G_out - Go_mask, party, io_server);
		COMM += 2 * G_out.rows() * G_out.cols() * sizeof(unsigned long int);
		W_grad = inner_product(in_x.transpose(), G_out, x_mask.transpose(), Go_mask, Wg_tp);
		truncation_mat(W_grad);
		bias_grad = G_out.colwise().sum();
		Mat rs = inner_product(G_out, W.transpose(), Go_mask, W_mask.transpose(), backward_tp);
		truncation_mat(rs);
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
		//cout << "done\n";
		return ans - Y;
	}

};
