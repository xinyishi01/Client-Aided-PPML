#ifndef NN_SERVER_HPP
#define NN_SERVER_HPP

#include "utils.hpp"
#include "Layer1_Server.hpp"
#include "Layer_Server.hpp"

class NN_Server{
	int party;
	
	Layer1_Server layer1;
	Layer_Server layer2;
	Layer_Server layer3;
	
	PRF prf_server, prf_client_same;
	
	PRF *prf_client;
	
	Mat train_data, train_label, X_batch, Y_batch;
	
	double accuracy_max, accuracy_now;
	
	double COMM, COMM_client;
	
	NetIO * io_server;
	NetIO * io_client;
	
public:

	
	NN_Server(int _party, NetIO *_io_server, NetIO *_io_client) {
		party = _party;
		
		io_server = _io_server;
		io_client = _io_client;
		
		if(party == ALICE) {
			prf_server.random_key();
			send_keys(&prf_server, io_server, 1);
		}
		else receive_keys(&prf_server, io_server, 1);
		
		layer1 = Layer1_Server(party, D, 128, "L1", learning_rate, &prf_server);
		layer2 = Layer_Server(party, 128, 128, "L2", learning_rate, &prf_server);
		layer3 = Layer_Server(party, 128, 10, "L3", learning_rate, &prf_server);
		
		layer1.set_io(io_server, io_client);
		layer2.set_io(io_server, io_client);
		layer3.set_io(io_server, io_client);
		
		train_data.resize(N, D);
		train_label.resize(N, K);
		X_batch.resize(B, D);
		Y_batch.resize(B, K);
		
		COMM = COMM_client = 0.;
	}
	
	void pre() {
		prf_client_same.random_key();
		prf_client = new PRF[B];
			
		COMM_client += send_keys(prf_client, io_client, B);
		COMM_client += send_keys(&prf_client_same, io_client, 1);
		
		layer1.pre(prf_client, &prf_client_same, X_batch);
		layer2.pre(prf_client, &prf_client_same);
		layer3.pre(prf_client, &prf_client_same);
	}
	
	Mat forward() {
		//cout << "l1...\n";
		Mat X1 = layer1.forward();
        if(DEBUG_FORWARD) {
        	Mat rec_X1 = reconstruct(X1, party, io_server);
        	if(party == ALICE) debug_mat(rec_X1, "X1.txt");
        }
        
		//cout << "l2...\n";
        Mat X2 = layer2.forward(X1);
        if(DEBUG_FORWARD) {
        	Mat rec_X2 = reconstruct(X2, party, io_server);
        	if(party == ALICE) debug_mat(rec_X2, "X2.txt");
        }
        
		//cout << "l3...\n";
		Mat X3 = layer3.forward(X2);
		if(DEBUG_FORWARD) {
        	Mat rec_X3 = reconstruct(X3, party, io_server);
        	if(party == ALICE) debug_mat(rec_X3, "X3.txt");
        }
        return X3;
	}
	
	void backward(Mat &X3, int it) {
		Mat G0 = layer3.softmax(X3, Y_batch);
		Mat G1 = layer3.backward(G0);
		Mat G2 = layer2.backward(G1);
		layer1.backward(G2);
	}
	
	void update() {
		layer1.update_W();
		layer2.update_W();
		layer3.update_W();
	}
	
	void evaluate(const Mat& test_data, const Mat& test_label, int it) {
			Mat W1 = reconstruct(layer1.W, party, io_server);
        	Mat bias1 = reconstruct(layer1.bias, party, io_server);
        	Mat W2 = reconstruct(layer2.W, party, io_server);
        	Mat bias2 = reconstruct(layer2.bias, party, io_server);
        	Mat W3 = reconstruct(layer3.W, party, io_server);
        	Mat bias3 = reconstruct(layer3.bias, party, io_server);
        	if(party == ALICE) {
        		
        		if(DEBUG_BACKWARD) {
        			debug_mat(W1, "W1.txt");
        			debug_mat(W2, "W2.txt");
        			debug_mat(W3, "W3.txt");
        			debug_mat(bias1, "b1.txt");
        			debug_mat(bias2, "b2.txt");
        			debug_mat(bias3, "b3.txt");
        		}
        		
		    	Mat X1 = plain_ReLU(test_data, W1, bias1);
				Mat X2 = plain_ReLU(X1, W2, bias2);
				Mat X3 = plain_ReLU(X2, W3, bias3);
        		
        		double tt = 0, hs = 0.001;
				for(int i = 0; i < X3.rows(); ++i) for (int j = 0; j < X3.cols(); ++j) {
					if((long int)X3(i, j) > 0) {
						tt += (long int)X3(i, j) / pow(2., L);
						++hs;
					}
				}
				
				//Mat X2 = plain_ReLU(test_data, W1, bias1);
				//Mat X3 = plain_ReLU(X2, W2, bias2);
			   	int ac = 0;
			   	for (int i = 0; i < testN; ++i) {
					long int maxvalue = X3(i, 0);
					int label = 0;
					for (int j = 1; j < K; ++j) {
						if((long int)X3(i, j) > maxvalue) {
							maxvalue = (long int)X3(i, j);
							label = j;
						}
					}
					if(label == test_label(i, 0)) ++ac;
				}
				accuracy_now = 100.0 * ac / testN;
				accuracy_max = max(accuracy_max, accuracy_now);
				cout << "#" << it <<": \t" << accuracy_now << "% \t" << accuracy_max << "%\n";
				cerr << "................................" << tt / hs << " " << (int)hs << endl;
			}
	}
	
	void online(const Mat& test_data, const Mat& test_label, vector<int>& perm) {
	
		cerr << "online ...\n";
		Time_Calculator mytime;
		double total_time = 0.;
		int t1 = mytime.start();
	
		prf_client = new PRF[N];
		COMM_client += send_keys(prf_client, io_client, N);
	
		COMM_client += merge_mat(train_data, prf_client, io_client, party);
		COMM_client += merge_mat(train_label, prf_client, io_client, party);
		
		delete [] prf_client;
		
		total_time += mytime.end(t1);
		
		double update_time = 0.;
		double pre_time = 0.;
		for (int it = 0, start = 0; it < IT; ++it, start += B) {
		
			if(it % 5 == 0) cout << "Iteration "<< it << endl;
			next_batch(X_batch,start,perm,train_data);
		    next_batch(Y_batch,start,perm,train_label);
		    int t1 = mytime.start();
			
			//cerr << "pre ...\n";
			int t2 = mytime.start();
			pre();
			pre_time += mytime.end(t2);
			
			//cerr << "forward ...\n";
			Mat X3 = forward();
			
			//cerr << "backward ...\n";
			backward(X3, it);
			
			t2 = mytime.start();
			//cerr << "update ...\n";
			update();
			update_time += mytime.end(t2);
			
			if(EVALUATE && it % 2 == 0) evaluate(test_data, test_label, it);
			
			delete[] prf_client;
			
			total_time += mytime.end(t1);
		}
		
		COMM += layer1.COMM + layer2.COMM + layer3.COMM;
		COMM_client += layer1.COMM_client + layer2.COMM_client + layer3.COMM_client;
		layer1.print();
		layer2.print();
		layer3.print();
		cout << "...............................\n";
		cout << "total time: " << total_time << endl;
		cout << "pre time: " << pre_time << endl;
		cout << "update time: " << update_time << endl;
		cout << "communication between servers:" << COMM / 1024 / 1024 << "MB\n";
		cout << "communication with clients:" << 2 * COMM_client / 1024 / 1024 << "MB\n";
	}
	
};

#endif
