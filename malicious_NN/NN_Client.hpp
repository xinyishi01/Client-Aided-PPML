#ifndef NN_CLIENT_HPP
#define NN_CLIENT_HPP

#include "utils.hpp"
#include "Layer1_Client.hpp"
#include "Layer_Client.hpp"

class NN_Client{

	Layer1_Client layer1;
	Layer_Client layer2;
	Layer_Client layer3;
	
	PRF *prf_alice, *prf_bob, prf_alice_same, prf_bob_same;
	
	Mat X_batch, Y_batch;
	
	NetIO *io_alice, *io_bob;
	
public:

	
	NN_Client(NetIO *_io_alice, NetIO *_io_bob) {
	
		io_alice = _io_alice;
		io_bob = _io_bob;
		
		layer1 = Layer1_Client(D, 128);
		layer2 = Layer_Client(128, 128);
		layer3 = Layer_Client(128, 10);
		
		layer1.set_io(io_alice, io_bob);
		layer2.set_io(io_alice, io_bob);
		layer3.set_io(io_alice, io_bob);
		
		layer1.set_id("L1");
		layer2.set_id("L2");
		layer3.set_id("L3");
		
		X_batch.resize(B, D);
		Y_batch.resize(B, K);
	}
	
	void pre() {
		prf_alice = new PRF[B];
		prf_bob = new PRF[B];
		receive_keys(prf_alice, io_alice, B);
		receive_keys(prf_bob, io_bob, B);
		receive_keys(&prf_alice_same, io_alice, 1);
		receive_keys(&prf_bob_same, io_bob, 1);
		
		layer1.pre(X_batch, prf_alice, prf_bob, &prf_alice_same, &prf_bob_same);
		layer2.pre(prf_alice, prf_bob, &prf_alice_same, &prf_bob_same);
		layer3.pre(prf_alice, prf_bob, &prf_alice_same, &prf_bob_same);
	}
	
	void forward() {
		layer1.forward();
		layer2.forward();
		layer3.forward();
	}
	
	void backward() {
		layer3.backward();
	}
	
	void update() {
		layer1.update();
		layer2.update();
		layer3.update();
	}
	
	void online(Mat& train_data, Mat& train_label, vector<int>&perm) {
		prf_alice = new PRF[N];
		prf_bob = new PRF[N];
		receive_keys(prf_alice, io_alice, N);
		receive_keys(prf_bob, io_bob, N);
		
		distribute_mat(train_data, prf_alice, prf_bob, io_alice, io_bob);
		distribute_mat(train_label, prf_alice, prf_bob, io_alice, io_bob);
		
		delete [] prf_alice;
		delete [] prf_bob;
		
		for (int it = 0, start = 0; it < IT; ++it, start += B) {
		
			//cerr << "batch ...\n";
			next_batch(X_batch,start,perm,train_data);
		    next_batch(Y_batch,start,perm,train_label);
		    
		    //cerr << "pre ...\n";
			pre();
			//cerr << "forward ...\n";
			forward();
			
			backward();
			
			update();
			
			delete[] prf_alice;
			delete[] prf_bob;
		}
		
	}

};

#endif
