#include "utils.hpp"

typedef uint8_t* vec_t;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

double total_time, GT;
double COMM, COMM_client;

PRF prf_client[B], prf_client_same, prf_server, own_prf;
My_Buffer dt;

vec_t tenc[P];

Time_Calculator mytime;

Mat reconstruct(Mat X0, int party) {
    Mat X1, X;
    X1 = X0;
    if (party == ALICE) {
        send_mat(X0, io_server);
        receive_mat(X1, io_server);
    }
    else {
        receive_mat(X1, io_server);
        send_mat(X0, io_server);
    }
    X = X0 + X1;
    return X;
}

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

        temp1 = temp3/(double)pow(2,L+8);
        
        temp2 = (long int)y(i);

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

void exchange_comparisons(Mat &comparisons, PRF *prf) {
	int sz = (B + 63) >> 6;
	Mat temp(sz, 2);
	prf->cal(temp);
	for (int i = 0; i < B; ++i) {
		int x = (i >> 6), y = (i & 63);
		for (int j = 0; j < 2; ++j) comparisons(i, j) = (temp(x, j) >> y) & 1;
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

void activation_function_server(Mat& d, Mat select_comparison, Mat &z1, Mat &z2, int party) {
    PRP* prp = new PRP();

    Mat d1(B, 1), d2(B, 1);
    uint64 one_half = 1 << (L + 7); //1 -> 1 <<(L + 8)
    for (int i = 0; i < B; i++) {
        // given d = d1 + d2 (shares, not the ones below) we want to find
        // d + 1/2 > 0 <=> d1 + 1/2 > -d2 and
        // d - 1/2 < 0 <=> -d1 > d2 - 1/2
        // (d + 1/2) * b1 * b2 + (1-b2)
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
    COMM_client += merge_mat(tau1_mat, prf_client, io_client, party);
	

    if (party == BOB) d1 = -1 * d1;
    
    Mat res1 = secure_inner_product_1(d1, z1, b1, 1, party) + tau1_mat;
    
    for (int i = 0; i < B; ++i) if(select_comparison(i, 0)) res1(i) = d1(i) - res1(i);

	COMM_client += merge_mat(b2, prf_client, io_client, party);
	COMM_client += merge_mat(tau2_mat, prf_client, io_client, party);

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

int main(int argc, char** argv) {

    //assert(N % B == 0);

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
	io_client->send_data(&hahaha1, sizeof(int));
	io_client->flush();
	
    if(party == ALICE) {
    	send_keys(&prf_server, 1, &dt);
    	dt.send(io_server);
    }
    else receive_keys(&prf_server, io_server, 1);
    

    cout<<"start training..."<<endl;
    
    int t1 = mytime.start();
    
    PRF *prf_client_temp = new PRF[N];
    COMM_client += send_keys(prf_client_temp, N, &dt);
    dt.send(io_client);
    
	
	COMM_client += merge_mat(train_data, prf_client_temp, io_client, party);
	COMM_client += merge_mat(train_label, prf_client_temp, io_client, party);
	
	delete [] prf_client_temp;
	
	total_time = mytime.end(t1);
	
	
    //start training
    for (int it = 0, start = 0; it < IT; ++it, start+=B) {
        
        	next_batch(x_batch,start,perm,train_data);
            next_batch(y_batch,start,perm,train_label);

            int t2 = mytime.start();
            
    		exchange_comparisons(comparisons_batch, &prf_server);
            
            //cout << "random keys..." << endl;
            for (int i = 0; i < B; ++i) prf_client[i].random_key();
            prf_client_same.random_key();
            
            //cout << "send keys..." << endl;
            COMM_client += send_keys(prf_client, B, &dt);
		    COMM_client += send_keys(&prf_client_same, 1, &dt);
		    dt.send(io_client);
		    
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
            activation_function_server(omega, comparisons_batch, z0_act_batch, z1_act_batch, party);
            //cout << "activation function done ...\n";
            omega = omega - y_batch * (uint64)(1<<(L+8));


            delta = secure_inner_product_0(omega, z1_batch, x_batch, party);
            
            for (int i = 0; i < B; ++i) delta += tau1_batch.row(i).transpose();

            for (int i=0; i<delta.rows(); i++) {
                delta(i) = (long int)delta(i)/ ((long int)(1<<23)*B);
            }

            W = W - delta;

            total_time += mytime.end(t2);
            

            if (EVALUATE && it % 5 == 4) {
            	NetIO * io = io_server;
                vector<long int> W_temp(W.cols()*W.rows());
                Mat W1(D,1);
				//cout << "EVALUATE...\n";
                if(party==ALICE) {

                    io->recv_data(&W_temp[0],sizeof(uint64)*W_temp.size());

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

                    io->send_data(&W_temp[0],sizeof(uint64)*W_temp.size());
                    io->flush();
                }
            }
        
    }

    cout<<"total time:"<<total_time<<"s!!"<<endl;
    cout<<"communication between servers:"<< COMM/1024/1024 << " MB" << endl;
    cout<<"communication with clients:"<<2.0 * COMM_client/1024/1024 << " MB" << endl;

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
