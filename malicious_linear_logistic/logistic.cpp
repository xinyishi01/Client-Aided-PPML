#include "utils.hpp"

typedef uint8_t* vec_t;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

double total_time, COMM, COMM_client;

PRF prf_client[B], prf_client_same, prf_server, own_prf;

vec_t tenc[P];
/*-------------------Malicious----------------*/

Mat cA0[cN + 10], cB0[cN + 10], cC0[cN + 10];
Mat /*cA1[cN + 10], */cB1[cN + 10], cC1[cN + 10];
int cNum0(0), cNum1(0);
Mat cX0[cM + 10], cY0[cM + 10], cZ0[cM + 10];
Mat /*cX1[cM + 10], */cY1[cM + 10], cZ1[cM + 10];
bool malicious_fl(0);

//Mat b1_, b2_, tau1_mat_, tau2_mat_;

/*-------------------Malicious----------------*/

/*-------------------Malicious-act----------------*/

uint64 dd[(dNN << 1) + 10], dr[(dNN << 1) + 10], db[(dNN << 1) + 10], dz[(dNN << 1) + 10];
int dNum(0);

/*-------------------Malicious-act----------------*/

Time_Calculator mytime;
void load_test_data(Mat& test_data, Mat& test_label) {
    ifstream infile2( "../data/mnist_test.csv" );
    int i=0;

    cout<<"load testing data.......\n";

    while(infile2) {

        string s;
        if (!getline(infile2,s))
            break;
        istringstream ss(s);
        int temp;
        char c;

        //read label
        ss>>temp;
        ss>>c;

        //if(temp == 0 || temp == 1){
        test_label(i) = (temp!=0);


        //read data (last entry 1)
        for(int j=0; j<D-1; j++) {
            ss>>test_data(i,j);
            ss>>c;
        }

        test_data(i,D-1) = 1;
        i++;
        //}

    }

    test_data.conservativeResize(i, D);
    test_label.conservativeResize(i,1);

    infile2.close();

    return;
}

double test_model(Mat& W0, Mat& W1, Mat& x, Mat& y) {
    Mat y_,W;
    double temp1;
    long int temp2,temp3;

    W = W0+W1;
    y_ = x*W;

    int count = 0;

    for(int i=0; i<y.rows(); i++) {
        temp3 = (long int)y_(i);
        //temp3 = (temp3<<4);

        //if(temp3>conv<long int>(p)/2)
        //	temp3 = temp3-conv<long int>(p);

        temp1 = temp3/(double)pow(2,L+8);
        //temp1 = conv<long int>(y_[i][0])/(double)pow(2,L);
        temp2 = (long int)y(i);

        //if(temp2>conv<long int>(p)/2)
        //	temp2 = temp2-conv<long int>(p);

        temp1 += 0.5; //activation function

        if(temp1>0.5 && temp2 == 1) {
            count++;
        } else if(temp1<0.5 && temp2 == 0) {
            count++;
        }
    }
    return count/(double)y.rows();


}

Mat secure_inner_product_D(Mat&x, Mat& z, Mat& y, int party) {

    Mat alpha0(D, 1), alpha1(D, 1);

    if(party == ALICE) alpha0 = x - z.transpose();
    else alpha1 = x - z.transpose();

    if(party == ALICE) {
        COMM += send_mat(alpha0, io_server);
        COMM += receive_mat(alpha1, io_server);

    } else {
        COMM += receive_mat(alpha0, io_server);
        COMM += send_mat(alpha1, io_server);
    }

    Mat result = y * (alpha0 + alpha1);

    return result;
}

Mat secure_inner_product_0(Mat&x, Mat& z, Mat& y, int party) {

    Mat alpha0(B, 1), alpha1(B, 1);

    if(party == ALICE) {
    	alpha0 = x - z;
    } else {
        alpha1 = x - z;
    }

    if(party == ALICE) {
        COMM += send_mat(alpha0, io_server);
        COMM += receive_mat(alpha1, io_server);

    } else {
        COMM += receive_mat(alpha0, io_server);
        COMM += send_mat(alpha1, io_server);
    }

    Mat result = (y.transpose()) * (alpha0 + alpha1);

    return result;
}


Mat secure_inner_product_1(Mat&x, Mat& z, Mat& y, int cols, int party) {

    Mat alpha0(B, cols), alpha1(B, cols);

    if(party == ALICE) {
        for (int col = 0; col < cols; col++) {
            alpha0.col(col) = x - z.col(col);
        }
    } else {
        for (int col = 0; col < cols; col++) {
            alpha1.col(col) = x - z.col(col);
        }
    }

    if(party == ALICE) {
        COMM += send_mat(alpha0, io_server);
        COMM += receive_mat(alpha1, io_server);

    } else {
        COMM += receive_mat(alpha0, io_server);
        COMM += send_mat(alpha1, io_server);
    }

    Mat result = (alpha0 + alpha1).cwiseProduct(y);

    return result;
}

void exchange_comparisons(int row, int col, Mat &comparisons, PRF *prf) {
	int sz = (row + 63) >> 6;
	Mat temp(sz, col);
	prf->cal(temp);
	for (int i = 0; i < row; ++i) {
		int x = (i >> 6), y = (i & 63);
		for (int j = 0; j < col; ++j)
            comparisons(i, j) = (temp(x, j) >> y) & 1;
	}
}

uint64 get_element(uint64 x, unsigned a, bool b, uint64 p, uint64 r) {
    uint64 modified;
    if ((unsigned) b == 0) {
        if (a == 0) {
            modified = (x + 1) << p;
        } else {
            modified = (1UL << (P - 1)) + r; // non-overlapping value
        }
    } else {
        if (a == 1) {
            modified = x << p;
        } else {
            modified = (3UL << (P - 2)) + r; // represented as 10 in binary so non-overlapping value
        }
    }
    return modified;
}

vec_t decompose_long(uint64 x, unsigned int length) {
    vec_t repr = new uint8_t[sizeof(uint64) * 2];
    for (size_t i = 0; i < sizeof(uint64) * 2; i++) {
        repr[i] = (uint8_t) x;
        x = x >> 8;
    }
    return repr;
}

void compute_encoding(bool deb, vec_t* vec, uint64 x, bool b) {
    unsigned a;
    uint64 modified;
    Mat r(P, 1);
    own_prf.cal(r);
    for (int i = 0; i < P; i++) {
        a = x & 1;
        modified = get_element(x, a, b, i, r(i, 0) >> 2);
        if (deb)
            cerr << i << ' ' << a << ' ' << b << ' ' << modified << ' ' << (r(i, 0) >> 2) << endl;
        vec[i] = decompose_long(modified, P - i + 1);
        x = x >> 1;
    }
}

void permute_encoding(vec_t* vec, PRP* prp) {
    for (int i = 0; i < P; i++) {
        prp->permute_data(vec[i], sizeof(uint64) * 2);
    }
}

void cal_encodings(bool deb, uint64 share, bool comp, int party, PRP* prp, vec_t encoding[]) {
    
    if(DEBUG_ACT) {
    	if(party == ALICE) {
    		bool comp_;
    		io_server->recv_data(&comp_, sizeof(bool));
    		if(comp_ != comp) assert(0);
    	}
    	else {
    		io_server->send_data(&comp, sizeof(bool));
    		io_server->flush();
    	}
    }
    
    if (party == ALICE) comp = !comp;

    compute_encoding(deb, &encoding[0], share, comp);
    permute_encoding(&encoding[0], prp);
}

static bool comp_(const int &x, const int &y) {
    for (int i = 0; i < 16; i++) {
        if (tenc[x][i] < tenc[y][i])
            return 1;
        if (tenc[x][i] > tenc[y][i])
            return 0;
    }
    return 0;
}

void malicious_act_store(uint64 temp_d, uint64 temp_z, uint64 temp_b, uint64 temp_tau, bool select_comp, int party) {
	dr[dNum] = temp_z;
	dd[dNum] = temp_d;
	
	if(select_comp == 0) {
		db[dNum] = temp_b;
		dz[dNum] = temp_tau;
	}
	else {
		db[dNum] = (party == ALICE) - temp_b;
		dz[dNum] = temp_z - temp_tau;
	}
    dNum++;
}

void activation_function_server(Mat &d, Mat &select_comparison, Mat &z1, Mat &z2, int party) {
    PRP* prp = new PRP();

    Mat d1(B, 1), d2(B, 1);
    unsigned long one_half = 1 << (L + 7); //1 -> 1 <<(L + 8)
    for (int i = 0; i < B; i++) {
        // given d = d1 + d2 (shares, not the ones below) we want to find
        // d + 1/2 > 0 <=> d1 + 1/2 > -d2 and
        // d - 1/2 < 0 <=> -d1 > d2 - 1/2
        if (party == ALICE) {
            d1(i) = d(i) + one_half;
            d2(i) = -d(i);
        }
        if (party == BOB) {
            d1(i) = -d(i);
            d2(i) = d(i) - one_half;
        }
    }
    
    vec_t encodings[2][B][P];
    for (int i = 0; i < B; i++) {
    	uint64 temp_x;
    	bool now, temp_comp;
    	
    	temp_x = d1(i);
    	if(select_comparison(i, 0) && party == BOB) temp_x = temp_x + 1;
    	now = ((long int)d1(i) < 0);
    	temp_comp = select_comparison(i, 0) ^ now;
    	if(now) temp_x = -temp_x;
    	cal_encodings(0, temp_x, temp_comp, party, prp, encodings[0][i]);
    	
    	temp_x = d2(i);
    	if(select_comparison(i, 1) && party == BOB) temp_x = temp_x + 1;
    	now = ((long int)d2(i) < 0);
    	temp_comp = select_comparison(i, 1) ^ now;
    	if(now) temp_x = -temp_x;
    	cal_encodings(0, temp_x, temp_comp, party, prp, encodings[1][i]);
    }
    
    uint8_t enc_t[2 * B * P * 16];
    int id[P];
    int cnt = 0;

    for (int i = 0; i < 2; i++)
    	for (int j = 0; j < B; j++) {
            for (int k = 0; k < P; k++) {
                id[k] = k;
                tenc[k] = encodings[i][j][k];
            }
            sort(id, id + P, comp_);
    		for (int k = 0; k < P; k++)
    			for (int l = 0; l < 16; l++)
    				enc_t[cnt++] = encodings[i][j][id[k]][l];
        }
        
    io_client->send_data(&enc_t[0], sizeof(uint8_t) * 2 * B * P * 16);
    io_client->flush();
    COMM_client += sizeof(uint8_t) * 2 * B * P * 16;
    
    for (int i = 0; i < B; ++i) for (int j = 0; j < P; ++j) {
    	delete [] encodings[0][i][j];
    	delete [] encodings[1][i][j];
    }
    
    Mat b1(B, 1), b2(B, 1), tau1_mat(B, 1), tau2_mat(B,1);
    
    COMM_client += merge_mat(b1, prf_client, io_client, party);
    COMM_client += merge_mat(b2, prf_client, io_client, party);
    COMM_client += merge_mat(tau1_mat, prf_client, io_client, party);
    COMM_client += merge_mat(tau2_mat, prf_client, io_client, party);
    
    /*-------------------Malicious-act----------------*/
    for (int i = 0; i < B; i++) {
    	malicious_act_store(d1(i), z1(i, 0), b1(i, 0), tau1_mat(i, 0), select_comparison(i, 0), party);
    	malicious_act_store(d2(i), z2(i, 0), b2(i, 0), tau2_mat(i, 0), select_comparison(i, 1), party);
    }
    /*-------------------Malicious-act----------------*/

    if (party == BOB) d1 = -1 * d1;
    
    Mat res1 = secure_inner_product_1(d1, z1, b1, 1, party) + tau1_mat;
    for (int i = 0; i < B; ++i) if(select_comparison(i, 0)) res1(i) = d1(i) - res1(i);

    Mat res = secure_inner_product_1(res1, z2, b2, 1, party) + tau2_mat;
    for (int i = 0; i < B; ++i) if(select_comparison(i, 1)) res(i) = res1(i) - res(i);
    
    for (int i = 0; i < B; ++i) {
        if(select_comparison(i, 1) == 0) { //not b2
            if(party == ALICE) d(i) = res(i) + (1ULL << (L+8)) * (- b2(i) + 1);
            else d(i) = res(i) - (1ULL << (L+8)) * b2(i);
        } else d(i) = res(i) + (1ULL << (L+8)) * b2(i);
    }
    delete prp;
}


/*-------------------Malicious----------------*/

void output_mat(Mat &x) {
    cerr << x.rows() << ' ' << x.cols() << endl;
    
    for (int i = 0; i < x.rows(); i++) {
        for (int j = 0; j < x.cols(); j++)
            cerr << x(i, j) << " ";
        cerr << endl;
    }
    cerr << endl;
}

bool Triple_Verif_With_Opening(Mat &X0, Mat &Y0, Mat &Z0, Mat &Y0_, Mat &Z0_, int party) {
    Mat X1, Y1, Z1, Y1_, Z1_;
    X1 = X0;
    Y1 = Y0;
    Z1 = Z0;
    Y1_ = Y0_;
    Z1_ = Z0_;
    if (party == ALICE) {
        COMM += receive_mat(X1, io_server);
        COMM += receive_mat(Y1, io_server);
        COMM += receive_mat(Z1, io_server);
        COMM += receive_mat(Y1_, io_server);
        COMM += receive_mat(Z1_, io_server);
    }
    else {
        COMM += send_mat(X0, io_server);
        COMM += send_mat(Y0, io_server);
        COMM += send_mat(Z0, io_server);
        COMM += send_mat(Y0_, io_server);
        COMM += send_mat(Z0_, io_server);
    }
    if (party == ALICE) {
        Mat Z = (X0 + X1) * (Y0 + Y1);
        Mat Z_ = sp_product(X0 + X1, Y0_ + Y1_);
    //    cerr << "Z:" << Z.rows() << ' ' << Z.cols() << endl;
    //    cerr << "Z_:" << Z_.rows() << ' ' << Z_.cols() << endl;
        return (Z == (Z0 + Z1)) && (Z_ == (Z0_ + Z1_));
    }
    return 1;
}


void Triple_Verif_Using_Another_Without_Opening(Mat &p, Mat &q0, Mat &q1, Mat &C0, Mat &C1, Mat &X, Mat &Y0, Mat &Z0, Mat &Y1, Mat &Z1, Mat &d0, Mat &d1, int party) {
    if (party == ALICE) {
        d0 = -C0 + Z0 + X * q0 + p * Y0 + p * q0;
        d1 = -C1 + Z1 + sp_product(X, q1) + sp_product(p, Y1) + sp_product(p, q1);
    }
    else {
        d0 = -C0 + Z0 + X * q0 + p * Y0;
        d1 = -C1 + Z1 + sp_product(X, q1) + sp_product(p, Y1);
    }
}

Mat scale_product(Mat X, uint64 c) {
    Mat ret = X;
    for (int i = 0; i < X.rows(); i++)
        for (int j = 0; j < X.cols(); j++)
            ret(i, j) *= c;
    return ret;
}

void my_shuffle(int *ans, int n) {
    Mat x(n, 1);
    prf_server.cal(x);
    for (int i = 0; i < n; ++i) ans[i] = i;
    
	for (int i = n - 1; i >= 0; --i) {
		unsigned long int now = x(i, 0) % (i+1);
		if(i != now) swap(ans[i], ans[now]);
	}
}

bool check_malicious(int party) {
        /*
        TODO:
            server generate x0, y0; store in cX0[i], cY0[i];
            x0: 1*d;  y0: d*1;
            randomly chosen client caculate x0*y0;
            client send z to server; store in cZ0[i];
            z0: 1*1;

            server generate x1, y1; store in cX1[i], cY1[i];
            x1: d*1;  y1: 1*1;
            randomly chosen client caculate x1*y1;
            client send z to server; store in cZ1[i];
            z1: d*1;
        */
    
    PRF *prf = new PRF;
    COMM_client += send_keys(prf, io_client, 1);
    
    for (int i = 0; i < cM; ++i) {
        cX0[i].resize(B, D); cY0[i].resize(D, 1);
        prf->cal(cX0[i]); prf->cal(cY0[i]);
        cY1[i].resize(1, B);
        prf->cal(cY1[i]);
    }

    Mat temp0(B, cM), temp1(cM, D);
    COMM_client += merge_mat_same(temp0, prf, io_client, party);
    for (int i = 0; i < cM; ++i) {
        cZ0[i] = temp0.col(i); // B*1
    }
    
    COMM_client += merge_mat_same(temp1, prf, io_client, party);
    for (int i = 0; i < cM; ++i) {
        cZ1[i] = temp1.row(i); // 1*d
    }
    
    delete prf;
    int id[cM];
    my_shuffle(id, cM);

    for (int i = ck * cN; i < cM; i++) {
        if (!Triple_Verif_With_Opening(cX0[id[i]], cY0[id[i]], cZ0[id[i]], cY1[id[i]], cZ1[id[i]], party))
            return false;
    }


    Mat sa0(B, 1), sa1(B, 1), sb0(1, D), sb1(1, D), d0, d1;
    uint64 c;
    Mat cm(cN * ck, 1);
    prf_server.cal(cm);

    int cnt(0);
    Mat p0[cN * ck], q00[cN * ck], q10[cN * ck];
    Mat p1[cN * ck], q01[cN * ck], q11[cN * ck];
    Mat p[cN * ck], q0[cN * ck], q1[cN * ck];
    for (int i = 0; i < cN; i++) {
        for (int j = 0; j < ck; j++) {
            p1[cnt] = p0[cnt] = cA0[i] - cX0[id[i * ck + j]];
            q01[cnt] = q00[cnt] = cB0[i] - cY0[id[i * ck + j]];
            q11[cnt] = q10[cnt] = cB1[i] - cY1[id[i * ck + j]];
            cnt++;
        }
    }
    if (party == ALICE) {
        for (int i = 0; i < cnt; i++) {
            COMM += send_mat(p0[i], io_server);
            COMM += send_mat(q00[i], io_server);
            COMM += send_mat(q10[i], io_server);
        }
        for (int i = 0; i < cnt; i++) {
            COMM += receive_mat(p1[i], io_server);
            COMM += receive_mat(q01[i], io_server);
            COMM += receive_mat(q11[i], io_server);
        }
    }
    else {
        for (int i = 0; i < cnt; i++) {
            COMM += receive_mat(p1[i], io_server);
            COMM += receive_mat(q01[i], io_server);
            COMM += receive_mat(q11[i], io_server);
        }
        for (int i = 0; i < cnt; i++) {
            COMM += send_mat(p0[i], io_server);
            COMM += send_mat(q00[i], io_server);
            COMM += send_mat(q10[i], io_server);
        }
    }
    for (int i = 0; i < cnt; i++) {
        p[i] = p0[i] + p1[i];
        q0[i] = q00[i] + q01[i];
        q1[i] = q10[i] + q11[i];
    }

    for (int i = 0; i < cN; i++) {
        for (int j = 0; j < ck; j++) {
            Triple_Verif_Using_Another_Without_Opening(p[i * ck + j], q0[i * ck + j], q1[i * ck + j], cC0[i], cC1[i],
                 cX0[id[i * ck + j]], cY0[id[i * ck + j]], cZ0[id[i * ck + j]], cY1[id[i * ck + j]], cZ1[id[i * ck + j]], d0, d1, party);
            c = cm(i * ck + j, 0);
            if (!i && !j) {
                sa0 = scale_product(d0, c);
                sb0 = scale_product(d1, c);
            }
            else {
                sa0 = sa0 + scale_product(d0, c);
                sb0 = sb0 + scale_product(d1, c);
            }
        }
    }

    if (party == ALICE) {
        COMM += receive_mat(sa1, io_server);
        COMM += receive_mat(sb1, io_server);
    }
    else {
        COMM += send_mat(sa0, io_server);
        COMM += send_mat(sb0, io_server);
    }

    if (party == ALICE) {
        return (-sa0 == sa1) && (-sb0 == sb1);
    }

    return true;
}

/*-------------------Malicious----------------*/


/*-------------------Malicious-act----------------*/

uint64 linear_com(Mat x) {
	Mat y(x.rows(), x.cols());
	prf_server.cal(y);
	return (x.cwiseProduct(y)).sum();
}

bool check_malicious_activation_function(int dN, int party) {
		PRP* prp = new PRP();
		int dM = dN * dk + dk + 1;
		
		Mat temp_comparisons(dM, 1), temp_b(dM, 1), temp_r(dM, 1), temp_z(dM, 1);
		exchange_comparisons(dM, 1, temp_comparisons, &prf_server);
		
		Mat mask(dM, 1);
		prf_server.cal(mask);
		int id[dM];
    	my_shuffle(id, dM);
    	vec_t* encodings[dM];
		for (int i = 0; i < dM; ++i){
			uint64 tempx;
			encodings[i] = new vec_t[P];
			if(id[i] >= dN * dk) tempx = mask(id[i], 0) + (party != ALICE) * (1 + rand() % 10);
			else tempx = dd[id[i] / dk] + mask(id[i], 0);
			
			if(temp_comparisons(id[i], 0) && party == BOB) tempx = tempx + 1;
			bool now = ((long int)tempx < 0);
			if(now) tempx = -tempx;
			bool temp_comp = temp_comparisons(id[i], 0) ^ now;
			cal_encodings(0, tempx, temp_comp, party, prp, encodings[i]);
		}
		
		int id2[P];
		uint8_t* enc_t = new uint8_t[dM * P * 16];
		for (int i = 0, cnt = 0; i < dM; i++) {
			for (int k = 0; k < P; k++) {
                id2[k] = k;
                tenc[k] = encodings[i][k];
            }
            sort(id2, id2 + P, comp_);
            
			for (int k = 0; k < P; k++) for (int l = 0; l < 16; l++)
				enc_t[cnt++] = encodings[i][id2[k]][l];
		}		
		for (int i = 0; i < dM; ++i) {
			for (int k = 0; k < P; ++k) delete[] encodings[i][k];
			delete[] encodings[i];
		}
		
		delete prp;
		
		io_client->send_data(&enc_t[0], sizeof(uint8_t) * dM * P * 16);
		COMM_client += sizeof(uint8_t) * dM * P * 16;
		io_client->flush();
		
		delete[] enc_t;
		
		prf_client_same.cal(temp_r);
		
		COMM_client += merge_mat_same(temp_b, &prf_client_same, io_client, party);
		COMM_client += merge_mat_same(temp_z, &prf_client_same, io_client, party);
		
		Mat rs(dM, 1), rs2(dM, 1), now_r(dN * dk, 1);
		
		for (int i = 0; i < dM; ++i) {
			if(temp_comparisons(id[i], 0)) { // flip back
				temp_b(i, 0) = (party == ALICE) - temp_b(i, 0);
				temp_z(i, 0) = temp_r(i, 0) - temp_z(i, 0);
			}
			if(id[i] >= dN * dk) {
				rs2(id[i], 0) = temp_b(i, 0);
				rs(id[i], 0) = temp_z(i, 0);
			}
			else {
				int o = id[i] / dk;
				rs2(id[i], 0) = db[o] - temp_b(i, 0);
				now_r(id[i], 0) = dr[o] - temp_r(i, 0);
				rs(id[i], 0) = temp_z(i, 0) - dz[o];
			}
		}
		
		Mat open_r = reconstruct(now_r, party, io_server);
		COMM += now_r.rows() * now_r.cols() * sizeof(unsigned long int);
		
		for (int i = 0; i < dM; ++i) if(id[i] < dN * dk) {
			rs(id[i], 0) += temp_b(i, 0) * open_r(id[i], 0);
		}
		
		uint64 ans = linear_com(rs) + linear_com(rs2), ans_;

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
		COMM += sizeof(uint64) * 2;
		
		return ans == -ans_;
		
}

/*-------------------Malicious-act----------------*/

int main(int argc, char** argv) {

    assert(N % B == 0);

    //setup connection
    int port, party;
    parse_party_and_port(argv, &party, &port);
    setup(party, io_server, io_client_alice, io_client_bob, io_client);

    total_time = 0.0;
    COMM = 0.0;
    COMM_client = 0.0;

    Mat test_data(testN,D), test_label(testN,1);

    if(party == ALICE && EVALUATE) {
        load_test_data(test_data, test_label);
    }

    cout<<"declare some matrices..."<<endl;

    Mat W(D,1);
    W.setZero();
    
    Mat x_batch(B,D), y_batch(B,1), z0_batch(1,D), z1_batch(B,1), tau0_batch(B,1), tau1_batch(B, D);
    Mat z0_act_batch(B, 1), z1_act_batch(B, 1);
    Mat comparisons_batch(B, 2);
    Mat train_data(N, D), train_label(N, 1);

    Mat update0(B, 1);
    Mat omega(B, 1);

    Mat update1(B, D);

    Mat delta(D, 1);
    
    vector<int> perm(N * Ep, 0);
    io_client->recv_data(&perm[0], sizeof(int) * perm.size());
	
	
	int hahaha1 = 666, hahaha2;
	io_client->send_data(&hahaha1, sizeof(int));
	io_client->flush();
	if(party == ALICE) {
		io_server->send_data(&hahaha1, sizeof(int));
		io_server->flush();
		io_server->recv_data(&hahaha2, sizeof(int));
	}
	else if(party == BOB) {
		io_server->recv_data(&hahaha2, sizeof(int));
		io_server->send_data(&hahaha1, sizeof(int));
		io_server->flush();
	}
	
    if(party == ALICE) COMM += send_keys(&prf_server, io_server, 1);
    else COMM += receive_keys(&prf_server, io_server, 1);

	
    cout<<"start training..."<<endl;
    
    int t1 = mytime.start();
    
	PRF *prf_client_temp = new PRF[N];
    COMM_client += send_keys(prf_client_temp, io_client, N);
    
	COMM_client += merge_mat(train_data, prf_client_temp, io_client, party);
	COMM_client += merge_mat(train_label, prf_client_temp, io_client, party);
	
	total_time += mytime.end(t1);
	
	delete [] prf_client_temp;
	
    //start training
    for (int round = 0, start = 0; round < Ep; round++) {

        for(int i=0; i<N/B; i++, start += B) {
        
			next_batch(x_batch,start,perm,train_data);
		    next_batch(y_batch,start,perm,train_label);
		    
            int t1 = mytime.start();

    		exchange_comparisons(B, 2, comparisons_batch, &prf_server);
            
            //cout << "random keys..." << endl;
            for (int i = 0; i < B; ++i) prf_client[i].random_key();
            prf_client_same.random_key();
            
            //cout << "send keys..." << endl;
            COMM_client += send_keys(prf_client, io_client, B);
		    COMM_client += send_keys(&prf_client_same, io_client, 1);
		    
            prf_client_same.cal(z0_batch);
            Mat temp(1, 1);
            for (int j = 0; j < B; ++j) {
            	prf_client[j].cal(temp);
            	z1_batch(j, 0) = temp(0, 0);
            }
		    
		    //cout << "receive triples..." << endl;
		    COMM_client += merge_mat(tau0_batch, prf_client, io_client, party);
		    COMM_client += merge_mat(tau1_batch, prf_client, io_client, party);
		    
		    prf_client_same.cal(z0_act_batch);
		    prf_client_same.cal(z1_act_batch);
		    
		    if(DEBUG_PRF) {
		    	if(party == ALICE) {
		    		Mat _tau1(B, D), tau1(B, D);
		    		Mat correct_tau1(B, D);
		    		receive_mat(_tau1, io_server);
		    		receive_mat(correct_tau1, io_client);
		    		tau1 = _tau1 + tau1_batch;
		    		for (int i = 0; i < B; ++i) for (int j = 0; j < D; ++j) 
						if(tau1(i, j) != correct_tau1(i, j)) {
							cerr << "prf wrong!\n";
							assert(0);
						}
		    	}
		    	else {
		    		send_mat(tau1_batch, io_server);
		    	}
		    }
		    
            
            update0 = secure_inner_product_D(W, z0_batch, x_batch, party);
            
            omega = update0 + tau0_batch;
            //cout << "activation function ...\n";

            /*-------------------Malicious----------------*/

            activation_function_server(omega, comparisons_batch, z0_act_batch, z1_act_batch, party);
            if (malicious_fl)
                break;

            /*-------------------Malicious----------------*/
            
            omega = omega - y_batch * (unsigned long int)(1<<(L+8));

            delta = secure_inner_product_0(omega, z1_batch, x_batch, party);
            for (int i = 0; i < B; ++i) delta += tau1_batch.row(i).transpose();

            for (int i=0; i<delta.rows(); i++) {
                delta(i) = (long int)delta(i)/ ((long int)(1<<23)*B);
            }

            W = W - delta;
            
            /*-------------------Malicious----------------*/

            cA0[cNum0] = x_batch; // B*d
            cB0[cNum0] = z0_batch.transpose(); // d*1
            cC0[cNum0] = tau0_batch; // B*1
            cNum0++;

            cB1[cNum1] = z1_batch.transpose(); // 1*B
            cC1[cNum1] = tau1_batch.colwise().sum(); // 1*d
            cNum1++;

            if (cNum0 == cN) {
            //    cout << "check_malicious\n";
                if (!check_malicious(party)) {
                    cout << "have found malicious\n";
                    malicious_fl = 1;
                    break;
                }
                cNum0 = cNum1 = 0;
            }

            /*-------------------Malicious----------------*/

            /*-------------------Malicious-act----------------*/

        //    cerr << "dNum: " << dNum << endl;
            if (dNum >= dNN) {
            //   	cerr << "check malicious activation function...\n";
               	
                if (!check_malicious_activation_function(dNum, party)) {
                    cout << "have found malicious while doing activation function\n";
                    malicious_fl = 1;
                    break;
                }
                dNum = 0;
            }

            /*-------------------Malicious-act----------------*/
            
            total_time += mytime.end(t1);

            if (EVALUATE && i % 5 == 4) {
            	NetIO * io = io_server;
                vector<long int> W_temp(W.cols()*W.rows());
                Mat W1(D,1);
				//cout << "EVALUATE...\n";
                if(party==ALICE) {

                    io->recv_data(&W_temp[0],sizeof(unsigned long int)*W_temp.size());

                    for(int j=0; j<W1.rows(); j++) {
                        for(int k=0; k<W1.cols(); k++)
                            W1(j,k) = W_temp[j*W1.cols()+k];
                    }

                    double res = test_model(W, W1, test_data, test_label);
                    cout<<res<<endl;
                } else {
                    for(int j=0; j<W.rows(); j++) {
                        for(int k=0; k<W.cols(); k++)
                            W_temp[j*W.cols()+k] = W(j,k);
                    }

                    io->send_data(&W_temp[0],sizeof(unsigned long int)*W_temp.size());
                    io->flush();
                }
            }
            /*-------------------Malicious----------------*/
            if (malicious_fl)
                break;
            /*-------------------Malicious----------------*/
        }
    }

    cout<<"total time:"<<total_time<<"s!!"<<endl;
    cout << "communication between servers: " << COMM/1024/1024 << " MB\n";
    cout << "communication with clients: " << 2.0*COMM_client/1024/1024 << " MB\n";

    if (party != EVE) {
        delete io_server;
    }
    if (party != BOB) {
        delete io_client_alice;
    }
    if (party != ALICE) {
        delete io_client_bob;
    }
}
