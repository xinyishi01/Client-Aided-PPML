#include "utils.hpp"

double total_time, online_time, offline_time, communication, computation, com_client;
double COMM, COMM_client, online_COMM_client;

Time_Calculator mytime;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

Mat reconstruct(Mat A, int party){
    vector<unsigned long int> A_temp(A.cols()*A.rows());
    
    if(party==ALICE){
        io_server->recv_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
        Mat A_(A.rows(),A.cols());
        for(int i=0;i<A_.rows();i++){
            for(int j=0;j<A_.cols();j++){
                A_(i,j) = A_temp[i*A_.cols()+j];
            }
        }
        
        Mat A_rec = A+A_;
        
        for(int i=0;i<A.rows();i++){
            for(int j=0;j<A.cols();j++)
                A_temp[i*A.cols()+j] = A(i,j);
        }
        io_server->send_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
        io_server->flush();

        COMM += (double)sizeof(unsigned long int)*A_temp.size();
        
        return A_rec;
    }
    else{
        for(int i=0;i<A.rows();i++){
            for(int j=0;j<A.cols();j++)
                A_temp[i*A.cols()+j] = A(i,j);
        }
        io_server->send_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
        io_server->flush();

        COMM += (double)sizeof(unsigned long int)*A_temp.size();
        
        io_server->recv_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
        Mat A_(A.rows(),A.cols());
        for(int i=0;i<A_.rows();i++){
            for(int j=0;j<A_.cols();j++){
                A_(i,j) = A_temp[i*A_.cols()+j];
            }
        }
        
        Mat A_rec = A+A_;
        
        return A_rec;
    }
}

void load_test_data(Mat& test_data, Mat& test_label){
	ifstream infile2( "../data/mnist_test.csv" );
	int count1=0,count2=0;
	int i=0;
	
	//cout<<"load testing data.......\n";
	
	while(infile2){
		
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
			for(int j=0;j<D-1;j++){
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

Mat spe_product(Mat a, Mat b) {
	Mat c(B, 1);
	for (int i = 0; i < B; i++) {
		c(i, 0) = a.row(i) * b.col(i);
	}
	return c;
}

Mat compute_delta0(Mat& W0,Mat& x0, Mat& y0, Mat& e, Mat& b0,Mat& c0, Mat& Wb1 ){
	Mat f,d;
	f = W0-b0+Wb1;
	d = x0*f+e*W0+c0-y0*(unsigned long int)(1<<(L+8));
	return d;
	
}

Mat compute_delta1(Mat& W1,Mat& x1, Mat& y1, Mat& e, Mat& b1,Mat& c1, Mat& Wb0 ){
	Mat f,d;
	f = W1-b1+Wb0;
	d = x1*f+e*(W1-f)+c1-y1*(unsigned long int)(1<<(L+8));
	return d;
	
}

void update_W0(Mat& W0,Mat& x0,Mat& e, Mat& b0,Mat& d0, Mat& c0, Mat& db1){
	Mat f, delta;
	
	f = d0-b0+db1;
	delta = x0*f+e*d0+c0;
	
	//reconstruct(delta,1);
	
	for(int i=0;i<delta.rows();i++){
		delta(i) = (long int)delta(i)/ ((long int)(1<<23)*B);
	}
	//reconstruct(delta,1);
	
	W0 = W0 - delta;
	
}

void update_W1(Mat& W1,Mat& x1,Mat& e, Mat& b1,Mat& d1, Mat& c1, Mat& db0){
	Mat f, delta;
	
	f = d1-b1+db0;
	
	delta = x1*f+e*(d1-f)+c1;
	
	//reconstruct(delta,2);
	
	for(int i=0;i<delta.rows();i++){
		delta(i) = (long int)delta(i)/ ((long int)(1<<23)*B);
	}
	
	//reconstruct(delta,2);
	
	W1 = W1 - delta;

}

double test_model(Mat& W0, Mat& W1, Mat& x, Mat& y){
	Mat y_,W;
	double temp1;
	long int temp2,temp3;
	
	W = W0+W1;
	y_ = x*W;
	
	int count = 0;
	
	for(int i=0;i<y.rows();i++){
		temp3 = (long int)y_(i);
		//temp3 = (temp3<<4);
		
		//if(temp3>conv<long int>(p)/2)
		//	temp3 = temp3-conv<long int>(p);
		
		temp1 = temp3/(double)pow(2,L+8);
		//temp1 = conv<long int>(y_[i][0])/(double)pow(2,L);
		temp2 = (long int)y(i);
		
		//if(temp2>conv<long int>(p)/2)
		//	temp2 = temp2-conv<long int>(p);
		
		if(temp1>0.5 && temp2 == 1){
			count++;
		}
		else if(temp1<0.5 && temp2 == 0){
			count++;
		}
	}	
	
	//cout << "rows = " << (int)y.rows() << endl;
	return count/(double)y.rows();
}

void sync_server() {
//	cerr << "start sync_server\n";
	int t1 = 666, t2;
	io_client->send_data(&t1, sizeof(int));
	io_client->flush();
	io_client->recv_data(&t2, sizeof(int));
//	cerr << "end sync_server\n";
}

void sync_between_servers(int party) {
	int syn1 = 666, syn2;
	if(party == ALICE) {
		io_server->send_data(&syn1, sizeof(int));
		io_server->flush();
		io_server->recv_data(&syn2, sizeof(int));
	}
	else if(party == BOB) {
		io_server->recv_data(&syn2, sizeof(int));
		io_server->send_data(&syn1, sizeof(int));
		io_server->flush();
	}
}

int main(int argc, char** argv){
	//setup connection
	int port, party;
	parse_party_and_port(argv, &party, &port);
	setup(party, io_server, io_client_alice, io_client_bob, io_client);

	vector<int> perm(N * Ep, 0);

	io_client->recv_data(&perm[0], sizeof(int) * perm.size());


	total_time = 0.0;
	offline_time = 0.0;
	online_time = 0.0;
    communication = 0.0;
    computation = 0.0;
    //com_client = 0.0;

    COMM = 0.0;
    COMM_client = 0.0;
    online_COMM_client = 0.0;
	double COMM_ = 0.;

	//--------receive_client-----------------
	
	Mat train_data(N,D), train_label(N,1);
    Mat a(N,D), xa(N, D);/*, b_1(D,N*Ep),c_1(B,IT),b_2(B,IT),c_2(D,N*Ep)*/;
    
	receive_mat(a, io_client, COMM_);
	
	sync_between_servers(party);
	sync_server();

	int t1 = mytime.start();
	
    receive_mat(train_data, io_client, online_COMM_client);
    receive_mat(train_label, io_client, online_COMM_client);
	
	online_time += mytime.end(t1);
	
	cout << online_time << endl;
	
	t1 = mytime.start();
	xa = reconstruct(train_data - a, party);
	
	online_time += mytime.end(t1);
	
	cout << online_time << endl;
		

	//--------receive_client-----------------

	Mat test_data(testN,D), test_label(testN,1);
	
	if(party == ALICE){
		load_test_data(test_data, test_label);
	}
	
	//cout<<"reading multiplication triples......\n";

	//cout<<(reconstruct(train_data-a,party)==xa)<<"!!!!\n";
	
	Mat W(D,1);
	W.setZero();
	
	Mat x_batch(B,D), tx_batch,txa_batch, y_batch(B,1), xa_batch(B,D);
	//Mat b_IT1(D,1), b_IT2(B,1),c_IT1(B,1), c_IT2(D, 1);
	Mat *b_IT1 = new Mat[IT], *b_IT2 = new Mat[IT], *c_IT1 = new Mat[IT], *c_IT2 = new Mat[IT];
	
	Mat W0b(D,1),W1b(D,1),d0_buf(B,1),d1_buf(B,1);
	
	total_time = online_time + offline_time;
	
	for (int i = 0; i < IT; ++i) {
		b_IT1[i].resize(D, 1);
		b_IT2[i].resize(B, 1);
		c_IT1[i].resize(B, 1);
		c_IT2[i].resize(D, 1);
		receive_mat(b_IT1[i], io_client, COMM_);
		receive_mat(b_IT2[i], io_client, COMM_);
		receive_mat(c_IT1[i], io_client, COMM_);
		receive_mat(c_IT2[i], io_client, COMM_);

	}
	
	int start = 0;

	sync_between_servers(party);
	cout << "start training...\n";
	
	//start training
	for (int i = 0; i < IT; i++){
	//	cerr << "i: " << i << endl;

		Mat d0, d1, W_rec, d;
		
		next_batch(x_batch, start, perm, train_data);
		next_batch(y_batch, start, perm, train_label);
		tx_batch = x_batch.transpose();
		
		next_batch(xa_batch,start,perm,xa);
		txa_batch = xa_batch.transpose();
	
		//sync_between_servers(party);
		int t3 = mytime.start();

		t1 = mytime.start();
		int t2 = mytime.start();

		start+= B;
		//send and receive data
		
	//	cerr << "0\n";
		
		if(party == ALICE){
			W0b = W - b_IT1[i];
		}
		else{
			W1b = W - b_IT1[i];
		}

		computation += mytime.end(t1);
		
		
	//	cerr << "1\n";

		if(party == ALICE){
			t1 = mytime.start();

			vector<unsigned long int> temp0(W0b.cols()*W0b.rows()),temp1(W1b.cols()*W1b.rows());

			for(int j=0;j<W0b.rows();j++){
				for(int k=0;k<W0b.cols();k++)
					temp0[j*W0b.cols()+k] = W0b(j,k);
			}

			computation += mytime.end(t1);

			t1 = mytime.start();
			
			io_server->send_data(&temp0[0],sizeof(unsigned long int)*temp0.size());
			io_server->flush();
			io_server->recv_data(&temp1[0],sizeof(unsigned long int)*temp1.size());

			COMM += (double)sizeof(unsigned long int) * temp1.size();

			communication += mytime.end(t1);

			t1 = mytime.start();

			for(int j=0;j<W1b.rows();j++){
				for(int k=0;k<W1b.cols();k++)
					W1b(j,k) = temp1[j*W1b.cols()+k];
			}

			computation += mytime.end(t1);
			
		}
		else{
			t1 = mytime.start();

			vector<unsigned long int> temp0(W0b.cols()*W0b.rows()),temp1(W1b.cols()*W1b.rows());
	
			for(int j=0;j<W1b.rows();j++){
				for(int k=0;k<W1b.cols();k++)
					temp1[j*W1b.cols()+k] = W1b(j,k);
			}

			computation += mytime.end(t1);

			t1 = mytime.start();
			
			io_server->recv_data(&temp0[0],sizeof(unsigned long int)*temp0.size());
			io_server->send_data(&temp1[0],sizeof(unsigned long int)*temp1.size());
			io_server->flush();

			COMM += (double)sizeof(unsigned long int) * temp0.size();

			communication += mytime.end(t1);

			t1 = mytime.start();
		
			for(int j=0;j<W0b.rows();j++){
				for(int k=0;k<W0b.cols();k++)
					W0b(j,k) = temp0[j*W0b.cols()+k];
			}

			computation += mytime.end(t1);
		
		}

	//	cerr << "1.5\n";
		
		//train step1

		t1 = mytime.start();
		
		if(party == ALICE){
			d = compute_delta0(W,x_batch, y_batch, xa_batch, b_IT1[i],c_IT1[i], W1b );
			d0_buf = d-b_IT2[i];
			//reconstruct(d0,party);
		}
		else{
			d = compute_delta1(W,x_batch, y_batch, xa_batch, b_IT1[i],c_IT1[i], W0b );
			d1_buf = d-b_IT2[i];
			//reconstruct(d1,party);
		}

		computation += mytime.end(t1);

	//	cout<< "!!! " << (reconstruct(x_batch,party)*reconstruct(W,party) - reconstruct(y_batch, party)*(unsigned long int)(1<<(L+8)) == reconstruct(d,party)) <<endl;
		
		//send and receive data

	//	cerr << "1.7\n";

		if(party == ALICE){
			t1 = mytime.start();

			vector<unsigned long int> d_temp0(d0_buf.cols()*d0_buf.rows()),d_temp1(d1_buf.cols()*d1_buf.rows());
	
			for(int j=0;j<d0_buf.rows();j++){
				for(int k=0;k<d0_buf.cols();k++)
					d_temp0[j*d0_buf.cols()+k] = d0_buf(j,k);
			}

			computation += mytime.end(t1);

			t1 = mytime.start();

			io_server->send_data(&d_temp0[0],sizeof(unsigned long int)*d_temp0.size());
			io_server->flush();
			io_server->recv_data(&d_temp1[0],sizeof(unsigned long int)*d_temp1.size());
			
			COMM += (double)sizeof(unsigned long int) * d_temp1.size();

			communication += mytime.end(t1);

			t1 = mytime.start();

			for(int j=0;j<d1_buf.rows();j++){
				for(int k=0;k<d1_buf.cols();k++)
					d1_buf(j,k) = d_temp1[j*d1_buf.cols()+k];
			}

			computation += mytime.end(t1);
		}
		else{
			t1 = mytime.start();

			vector<unsigned long int> d_temp0(d0_buf.cols()*d0_buf.rows()),d_temp1(d1_buf.cols()*d1_buf.rows());
	
			for(int j=0;j<d1_buf.rows();j++){
				for(int k=0;k<d1_buf.cols();k++)
					d_temp1[j*d1_buf.cols()+k] = d1_buf(j,k);
			}

			computation += mytime.end(t1);

			t1 = mytime.start();
			
			io_server->recv_data(&d_temp0[0],sizeof(unsigned long int)*d_temp0.size());
			io_server->send_data(&d_temp1[0],sizeof(unsigned long int)*d_temp1.size());
			io_server->flush();
			
			COMM += (double)sizeof(unsigned long int) * d_temp0.size();

			communication += mytime.end(t1);

			t1 = mytime.start();

			for(int j=0;j<d0_buf.rows();j++){
				for(int k=0;k<d0_buf.cols();k++)
					d0_buf(j,k) = d_temp0[j*d0_buf.cols()+k];
			}

			computation += mytime.end(t1);
		
		}

	//	cerr << "2\n";
		
		//train step2

		t1 = mytime.start();
		
		if(party == ALICE){
			update_W0(W,tx_batch,txa_batch, b_IT2[i],d, c_IT2[i], d1_buf); 
		}
		else{
			update_W1(W,tx_batch, txa_batch, b_IT2[i],d, c_IT2[i], d0_buf);
		}

	//	cerr << "3\n";

		computation += mytime.end(t1);

		online_time += mytime.end(t2);

		total_time += mytime.end(t3);
		
		//reconstruct(W,party);
		
		if(EVALUATE){
			vector<long int> W_temp(W.cols()*W.rows());
			int tempi = 10;
			Mat W1(D,1);
			
			if(party==ALICE){
			
				io_server->send_data(&tempi,4);
				io_server->flush();
				io_server->recv_data(&W_temp[0],sizeof(unsigned long int)*W_temp.size());
				
				
				for(int j=0;j<W1.rows();j++){
					for(int k=0;k<W1.cols();k++)
						W1(j,k) = W_temp[j*W1.cols()+k];
				}
				
				
				double res = test_model(W,W1,test_data, test_label);
				cout<<res<<endl;
			}
			else{
				for(int j=0;j<W.rows();j++){
					for(int k=0;k<W.cols();k++)
						W_temp[j*W.cols()+k] = W(j,k);
				}
				
				io_server->recv_data(&tempi,4);
				io_server->send_data(&W_temp[0],sizeof(unsigned long int)*W_temp.size());
				io_server->flush();
				
			}
			
		}
		
	}
	
	cout << "total time:" << total_time << "s" << endl;
	cout << "online time:" << online_time << "s" << endl;
	cout << "offline_time:" << offline_time << "s" << endl;
	cout << "comm between servers(online): " << (COMM / 512 / 1024) << "MB" << endl;
	cout << "offline comm with clients: " << (COMM_client / 512 / 1024) << "MB" << endl;
	cout << "online comm with clients: " << (online_COMM_client / 512 / 1024) << "MB" << endl;
	cout << "\n";
	
	delete[] b_IT1;
	delete[] b_IT2;
	delete[] c_IT1;
	delete[] c_IT2;
	
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
