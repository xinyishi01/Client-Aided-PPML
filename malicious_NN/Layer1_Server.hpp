#include "utils.hpp"
#include "Layer.hpp"


class Layer1_Server : public Layer {
public:
	
	Mat z, tau, backward_z, backward_tau;

/*--------------------------Malicious--------------------------*/
	Mat oX[cN], oZ[cN], oTAU[cN], obZ[cN], obTAU[cN];
/*--------------------------Malicious--------------------------*/
	
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
		COMM = COMM_client = 0.;
		
		/*--------------------------Malicious--------------------------*/
		oNum = 0;
		dNum = 0;
		/*--------------------------Malicious--------------------------*/
	}


/*--------------------------Malicious--------------------------*/

bool check_malicious() {
    PRF *prf = new PRF;
    COMM_client += send_keys(prf, io_client, 1);
    
    Mat *nwX = new Mat[cM];
    Mat *nwZ = new Mat[cM];
    Mat *nwTAU = new Mat[cM];
    Mat *nwbZ = new Mat[cM];
    Mat *nwbTAU = new Mat[cM];
    
    for (int i = 0; i < cM; ++i) {
        nwX[i].resize(B, D_in); 
        nwZ[i].resize(D_in, D_out);
        nwTAU[i].resize(B, D_out);
        nwbZ[i].resize(B, D_out);
        nwbTAU[i].resize(D_in, D_out);
        prf->cal(nwX[i]); 
        prf->cal(nwZ[i]);
        prf->cal(nwbZ[i]);
    }

	Mat temp0(cM, B * D_out), temp1(cM, D_in * D_out);
	COMM_client += merge_mat_same(temp0, prf, io_client, party);
	COMM_client += merge_mat_same(temp1, prf, io_client, party);
	
	for (int i = 0; i < cM; ++i) {
		unfold_mat(nwTAU[i], temp0.row(i));
		unfold_mat(nwbTAU[i], temp1.row(i));
	}
	
	delete prf;
    
    int id[cM];
    my_shuffle(id, cM);
    bool fl = 1;
    
    uint64 ans = 0;
	
	Mat openX(nwX[0].rows(), nwX[0].cols());
	Mat openZ(nwZ[0].rows(), nwZ[0].cols());
	Mat openbZ(nwbZ[0].rows(), nwbZ[0].cols());
    for (int i = ck * cN; i < cM && fl; i++) {
        int o = id[i];/*
        Mat openX = reconstruct(nwX[o], party, io_server);
        Mat openZ = reconstruct(nwZ[o], party, io_server);
        Mat openTAU = reconstruct(nwTAU[o], party, io_server);
        Mat openbZ = reconstruct(nwbZ[o], party, io_server);
        Mat openbTAU = reconstruct(nwbTAU[o], party, io_server);
        if((openX * openZ) != openTAU) fl = 0;
        if((openX.transpose() * openbZ) != openbTAU) fl = 0;*/
        if(party == ALICE) {
        	COMM += receive_mat(openX, io_server);
        	COMM += receive_mat(openZ, io_server);
        	COMM += receive_mat(openbZ, io_server);
        	openX = openX + nwX[o];
        	openZ = openZ + nwZ[o];
        	openbZ = openbZ + nwbZ[o];
        	ans += linear_com(openX * openZ - nwTAU[o]);
        	ans += linear_com(openX.transpose() * openbZ - nwbTAU[o]);
        }
        else {
        	COMM += send_mat(nwX[o], io_server);
        	COMM += send_mat(nwZ[o], io_server);
        	COMM += send_mat(nwbZ[o], io_server);
        	ans += linear_com(-nwTAU[o]);
        	ans += linear_com(-nwbTAU[o]);
        }
        
    }
	
	//if(fl) {
	Mat *tempe = new Mat[cN * ck];
	Mat *tempf = new Mat[cN * ck];
	Mat *tempbf = new Mat[cN * ck];
	Mat *tempe_ = new Mat[cN * ck];
	Mat *tempf_ = new Mat[cN * ck];
	Mat *tempbf_ = new Mat[cN * ck];
	for (int i = 0, cnt = 0; i < cN; i++) {
	    for (int j = 0; j < ck; j++, ++cnt) {
	    	int o = id[i * ck + j];
	    	tempe[cnt] = oX[i] - nwX[o];
	    	tempf[cnt] = oZ[i] - nwZ[o];
	    	tempbf[cnt] = obZ[i] - nwbZ[o];
	    	tempe_[cnt].resize(tempe[cnt].rows(), tempe[cnt].cols());
	    	tempf_[cnt].resize(tempf[cnt].rows(), tempf[cnt].cols());
	    	tempbf_[cnt].resize(tempbf[cnt].rows(), tempbf[cnt].cols());
	   }
	}
	   
	if(party == ALICE) {
		for (int i = 0; i < cN * ck; ++i) {
			COMM += send_mat(tempe[i], io_server);
			COMM += send_mat(tempf[i], io_server);
			COMM += send_mat(tempbf[i], io_server);
		}
		for (int i = 0; i < cN * ck; ++i) {
			COMM += receive_mat(tempe_[i], io_server);
			COMM += receive_mat(tempf_[i], io_server);
			COMM += receive_mat(tempbf_[i], io_server);
		}
	}
	else {
		for (int i = 0; i < cN * ck; ++i) {
			COMM += receive_mat(tempe_[i], io_server);
			COMM += receive_mat(tempf_[i], io_server);
			COMM += receive_mat(tempbf_[i], io_server);
		}
		for (int i = 0; i < cN * ck; ++i) {
			COMM += send_mat(tempe[i], io_server);
			COMM += send_mat(tempf[i], io_server);
			COMM += send_mat(tempbf[i], io_server);
		}
	}
	
	for (int i = 0; i < cN * ck; ++i) {
		tempe[i] = tempe[i] + tempe_[i];
		tempf[i] = tempf[i] + tempf_[i];
		tempbf[i] = tempbf[i] + tempbf_[i];
	}
	
		Mat rs(cN * ck, 2);
		for (int i = 0, cnt = 0; i < cN; i++) {
		    for (int j = 0; j < ck; j++, ++cnt) {
		    	int o = id[i * ck + j];/*
		    	Mat tempe = reconstruct(oX[i] - nwX[o], party, io_server);
		    	Mat tempf = reconstruct(oZ[i] - nwZ[o], party, io_server);
		    	Mat tempbf = reconstruct(obZ[i] - nwbZ[o], party, io_server);*/
        		
		        rs(i * ck + j, 0) = linear_com(inner_product(oX[i], oZ[i], tempe[cnt], tempf[cnt], nwTAU[o]) - oTAU[i]);
		        rs(i * ck + j, 1) = linear_com(inner_product(oX[i].transpose(), obZ[i], tempe[cnt].transpose(), tempbf[cnt], nwbTAU[o]) - obTAU[i]);
		    }
		}
		ans += linear_com(rs);
	delete []tempe;
	delete []tempf;
	delete []tempbf;
	delete []tempe_;
	delete []tempf_;
	delete []tempbf_;
	
		uint64 ans_;
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
	
	COMM += 2 * sizeof(unsigned long int);
	
    delete[] nwX;
    delete[] nwZ;
    delete[] nwTAU;
    delete[] nwbZ;
    delete[] nwbTAU;
    
    return fl;
}
		
/*--------------------------Malicious--------------------------*/	


	void pre(PRF *_prf_client, PRF *_prf_client_same, Mat _in_x) {
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
		
		/*--------------------------Malicious--------------------------*/
		oX[oNum] = in_x; //B, D_in
		oZ[oNum] = z; //D_in, D_out
		oTAU[oNum] = tau; //B, D_out
		obZ[oNum] = backward_z;//B, D_out
		obTAU[oNum] = backward_tau;//D_in, D_out
		
		++oNum;
		
		if(oNum == cN) {
			//cout << id << " is checking malicious...\n";
			if(!check_malicious()) {
				cout << "have found malicious!\n";
				assert(0);
			}
			oNum = 0;
			//cout << "OK.\n";
		}
		/*--------------------------Malicious--------------------------*/
		
	}
	Mat forward() {
		
		Mat inner = inner_product_input(z, in_x, tau);
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
	
	void backward(Mat& x) {
		++it;
		Mat G_out(B, D_out);
		G_out = cwise_product(x, z2, b1, party) + tau2;
		for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) if(select_comparison(i, j)) 
			G_out(i, j) = x(i, j) - G_out(i, j);
		if(DEBUG_BACKWARD) {
			Mat rec_G_out = reconstruct(G_out, party, io_server);
			if(party == ALICE) debug_mat(rec_G_out, id + "_G_out.txt");
		}
		
		Mat alpha = reconstruct(G_out - backward_z, party, io_server);
		COMM += 2 * G_out.rows() * G_out.cols() * sizeof(unsigned long int);
		
		W_grad = (in_x.transpose() * alpha) + backward_tau;
		truncation_mat(W_grad);
		
		bias_grad = G_out.colwise().sum();
	}

};
