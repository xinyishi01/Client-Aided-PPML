#include "utils.hpp"

double total_time, online_time, offline_time, communication, computation, com_client;
double COMM, COMM_client, COMM_, online_COMM_client;

Time_Calculator mytime;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

SHOTIterated* ot;
SHOTIterated* ot2;

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
	
	int i=0;
	
	cout<<"load testing data.......\n";
	
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

Mat compute_inner0(Mat& W0,Mat& x0, Mat& e, Mat& b0,Mat& c0, Mat& Wb1 ){
	Mat f,d;
	f = Wb1 - b0;
	for (int i = 0; i < B; i++)
		f.col(i) = f.col(i) + W0; 
	d = spe_product(x0, f)+e*W0+c0;
	return d;
	
}

Mat compute_inner1(Mat& W1,Mat& x1, Mat& e, Mat& b1,Mat& c1, Mat& Wb0 ){
	Mat f, f_, d;
	f = f_ = Wb0 - b1;
	for (int i = 0; i < B; i++) {
		f.col(i) = f.col(i) + W1;
		f_.col(i) = W1 - f.col(i);
	}
	d = spe_product(x1, f)+ spe_product(e, f_)+c1;
	return d;
	
}

Mat logistic(Mat& inner, Mat& y, int party){
	Batcher batcher1, batcher2;
	
	if(party==ALICE){
		for(int i = 0; i < inner.rows(); ++i) {
			batcher1.add<Integer>(P, inner(i));
			batcher2.add<Integer>(P, (unsigned long)0);
		}
	}
	else{
		for(int i = 0; i < inner.rows(); ++i) {
			batcher1.add<Integer>(P, (unsigned long)0);
			batcher2.add<Integer>(P, inner(i));	
		}
	}
	batcher1.make_semi_honest(ALICE);
	batcher2.make_semi_honest(BOB);

	COMM += (double)2 * inner.rows() * LAM * (2 * P - 1) / 8 / 2;
	
	Bit* a = new Bit[inner.rows()];
	Bit* b = new Bit[inner.rows()];
	Bit* c = new Bit[inner.rows()];
	Bit* d = new Bit[inner.rows()];
	bool* lsb_c = new bool[inner.rows()];
	bool* lsb_d = new bool[inner.rows()];
	
	//Integer temp2(L+1,(int)pow(2,L),BOB);
	
	for(int i = 0; i < inner.rows(); ++i){
		//Integer temp = batcher1.next<Integer>() + batcher2.next<Integer>();
		Integer temp = batcher2.next<Integer>();
		Integer temp2 = batcher1.next<Integer>()+temp;
		b[i] = (temp2)[P-1];
		a[i] = (temp2 - Integer(P, 1ULL<<L, 0))[P-1] ;
		//a[i] = (batcher1.next<Integer>() + temp)[P-1];
		
		c[i] = a[i]&(!(b[i]));
		d[i] = !(a[i]);
		lsb_c[i] = getLSB(c[i]);
		lsb_d[i] = getLSB(d[i]);
	}
	
	bool* select;
	select = new bool[inner.rows()*2];
	copy(lsb_c,lsb_c+inner.rows(),select);
	copy(lsb_d,lsb_d+inner.rows(),select+inner.rows());
	
	block* m0,*m1,*r;
	
	
	if(party==ALICE){
		m0 = new block[inner.rows()*2];
		m1 = new block[inner.rows()*2];
		r = new block[inner.rows()];
	}
	else{
		m0 = new block[inner.rows()];
		m1 = new block[inner.rows()];
		r = new block[inner.rows()*2];
	}
	
	//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 1"<<endl;
	
	uint64 temp_m0, temp_m1;
	Mat r1(inner.rows(), inner.cols()),r2,R(inner.rows(), inner.cols()),R2;
	
	for(int i=0;i<r1.rows();i++){
		r1(i) = randomlong();
	}
	
	
	if(party==ALICE){
		r2.resize(inner.rows(), inner.cols());
		
		for(int i=0;i<r2.rows();i++){
			r2(i) = randomlong();
		}
	}
	else{
		R2.resize(inner.rows(), inner.cols());
	}	
	
	for(int i=0;i<inner.rows();++i){
		temp_m0 = (int)select[i]*inner(i)+r1(i);
		temp_m1 = (int)(!select[i])*inner(i)+r1(i);
		m0[i] = makeBlock((uint64)0,temp_m0);
		m1[i] = makeBlock((uint64)0,temp_m1);
	}
	
	if(party == ALICE){
		for(int i=inner.rows();i<inner.rows()*2;++i){
			temp_m0 = (int)select[i]*(uint64)(1<<L)+r2(i-inner.rows());
			temp_m1 = (int)(!select[i])*(uint64)(1<<L)+r2(i-inner.rows());
			m0[i] = makeBlock((uint64)(0),temp_m0);
			m1[i] = makeBlock((uint64)(0),temp_m1);
		}
	}
	

	Mat y_;
	
	if(party == ALICE){
	
		ot->send(m0,m1,inner.rows()*2);

		COMM += (double)inner.rows() * 2.0 * LAM * 2 / 8;
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 4"<<endl;
		
		ot2->recv(r,select,inner.rows());

		COMM += (double)inner.rows() * 2.0 * LAM / 8;
		
		for(int i=0;i<inner.rows();i++){	
			R(i) = ((unsigned long *)(&r[i]))[0];
		}
		
		y_ = R-r1-r2-y*(uint64)(1<<L);
		
	}
	else{
		ot->recv(r,select,inner.rows()*2);

		COMM += (double)inner.rows() * 2.0 * LAM * 2 / 8;
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 4"<<endl;
		
		ot2->send(m0,m1,inner.rows());

		COMM += (double)inner.rows() * 2.0 * LAM / 8;

		for(int i=0;i<inner.rows();i++){	
			R(i) = ((unsigned long *)(&r[i]))[0];
			R2(i) = ((unsigned long *)(&r[i+inner.rows()]))[0];
		}
		
		y_ = R+R2-r1-y*(uint64)(1<<L);
	}
	
	
	//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 3"<<endl;
	
	return y_;
	
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
		
		temp1 = temp3/(double)pow(2,L);
		temp2 = (long int)y(i);
		
		temp1 += 0.5; //active
		
		if(temp1>0.5 && temp2 == 1){
			count++;
		}
		else if(temp1<0.5 && temp2 == 0){
			count++;
		}
	}
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
	
	setup_semi_honest(io_server, party);
	ot = new SHOTIterated(io_server, party==ALICE);
	ot2 = new SHOTIterated(io_server, party!=ALICE);
	
	
	vector<int> perm(N * Ep, 0);

	io_client->recv_data(&perm[0], sizeof(int) * perm.size());

	sync_between_servers(party);
	sync_server();

	total_time = 0.0;
	offline_time = 0.0;
	online_time = 0.0;
    communication = 0.0;
    computation = 0.0;
    com_client = 0.0;
    

    COMM = 0.0;
    COMM_client = 0.0;
    online_COMM_client = 0.0;
    COMM_ = (double)sizeof(int) * perm.size();

	//--------receive_client-----------------
	
	Mat train_data(N,D), train_label(N,1);
    Mat a(N,D), xa(N, D);/*, b_1(D,N*Ep),c_1(B,IT),b_2(B,IT),c_2(D,N*Ep)*/;

    int t1 = mytime.start();

	receive_mat(a, io_client, COMM_client);

	offline_time += mytime.end(t1);

	t1 = mytime.start();
	
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
	Mat b_IT1(D,B), b_IT2(B,1),c_IT1(B,1), c_IT2_(D,B), c_IT2(D, 1);
	
	Mat W0b(D,B),W1b(D,B),d0_buf(B,1),d1_buf(B,1);

	total_time = online_time + offline_time;
	
	int start = 0;

	sync_between_servers(party);
	sync_server();
	COMM_ += (double)sizeof(int) * 3;
	
	//start training
	for(int i=0;i<IT;i++){
	//	cerr << "it: " << i << endl;
	
		Mat d0,d1,W_rec, inner;
		
		next_batch(x_batch,start,perm,train_data);
		next_batch(y_batch,start,perm,train_label);
		tx_batch = x_batch.transpose();
		
		next_batch(xa_batch,start,perm,xa);
		txa_batch = xa_batch.transpose();
		
		sync_between_servers(party);
		sync_server();
		COMM_ += (double)sizeof(int) * 6;

		int t3 = mytime.start();
		t1 = mytime.start();

		receive_mat(b_IT1, io_client, COMM_client);
		receive_mat(b_IT2, io_client, COMM_client);
		receive_mat(c_IT1, io_client, COMM_client);
		receive_mat(c_IT2_, io_client, COMM_client);

		offline_time += mytime.end(t1);

	//	sync_between_servers(party);

	//	t1 = mytime.start();
		int t2 = mytime.start();
		
		c_IT2 = c_IT2_.col(0);
		for (int i = 1; i < B; i++) {
			c_IT2 = c_IT2 + c_IT2_.col(i);
		}
		
		//Mat a_batch(B,D);
		//next_batch(a_batch,start,perm,a);
		//cout<<(reconstruct(a_batch,party)*reconstruct(b_IT1,party) == reconstruct(c_IT1,party))<<endl;
		
		start+= B;
		
		//send and receive data		
		if(party == ALICE){
			for (int i = 0; i < B; i++)
				W0b.col(i) = W;
			W0b = W0b - b_IT1;
		}
		else{
			for (int i = 0; i < B; i++)
				W1b.col(i) = W;
			W1b = W1b - b_IT1;
		}
		
		vector<unsigned long int> temp0(W0b.cols()*W0b.rows()),temp1(W1b.cols()*W1b.rows());
	
		if(party == ALICE){
	
			for(int j=0;j<W0b.rows();j++){
				for(int k=0;k<W0b.cols();k++)
					temp0[j*W0b.cols()+k] = W0b(j,k);
			}
			
			io_server->send_data(&temp0[0],sizeof(unsigned long int)*temp0.size());
			io_server->flush();
			io_server->recv_data(&temp1[0],sizeof(unsigned long int)*temp1.size());
			
			COMM += (double)sizeof(unsigned long int) * temp1.size();

			for(int j=0;j<W1b.rows();j++){
				for(int k=0;k<W1b.cols();k++)
					W1b(j,k) = temp1[j*W1b.cols()+k];
			}
			
		}
		else{
			for(int j=0;j<W1b.rows();j++){
				for(int k=0;k<W1b.cols();k++)
					temp1[j*W1b.cols()+k] = W1b(j,k);
			}
			
			io_server->recv_data(&temp0[0],sizeof(unsigned long int)*temp0.size());
			io_server->send_data(&temp1[0],sizeof(unsigned long int)*temp1.size());
			io_server->flush();

			COMM += (double)sizeof(unsigned long int) * temp0.size();
		
			for(int j=0;j<W0b.rows();j++){
				for(int k=0;k<W0b.cols();k++)
					W0b(j,k) = temp0[j*W0b.cols()+k];
			}
		
		}	
		
		//train step1
		
		if(party == ALICE){
			inner = compute_inner0(W,x_batch, xa_batch, b_IT1,c_IT1, W1b );
			d0 = logistic(inner, y_batch, party);
			d0_buf = d0-b_IT2;
			//reconstruct(d0,party);
		}
		else{
			inner = compute_inner1(W,x_batch, xa_batch, b_IT1,c_IT1, W0b );
			d1 = logistic(inner, y_batch, party);
			d1_buf = d1-b_IT2;
			//reconstruct(d1,party);
		}
		
		//send and receive data
		vector<unsigned long int> d_temp0(d0_buf.cols()*d0_buf.rows()),d_temp1(d1_buf.cols()*d1_buf.rows());
		
		if(party == ALICE){
			for(int j=0;j<d0_buf.rows();j++){
				for(int k=0;k<d0_buf.cols();k++)
					d_temp0[j*d0_buf.cols()+k] = d0_buf(j,k);
			}

			io_server->send_data(&d_temp0[0],sizeof(unsigned long int)*d_temp0.size());
			io_server->flush();
			io_server->recv_data(&d_temp1[0],sizeof(unsigned long int)*d_temp1.size());
			
			COMM += (double)sizeof(unsigned long int) * d_temp1.size();

			for(int j=0;j<d1_buf.rows();j++){
				for(int k=0;k<d1_buf.cols();k++)
					d1_buf(j,k) = d_temp1[j*d1_buf.cols()+k];
			}
		}
		else{
			for(int j=0;j<d1_buf.rows();j++){
				for(int k=0;k<d1_buf.cols();k++)
					d_temp1[j*d1_buf.cols()+k] = d1_buf(j,k);
			}
			
			io_server->recv_data(&d_temp0[0],sizeof(unsigned long int)*d_temp0.size());
			io_server->send_data(&d_temp1[0],sizeof(unsigned long int)*d_temp1.size());
			io_server->flush();

			COMM += (double)sizeof(unsigned long int) * d_temp0.size();

			for(int j=0;j<d0_buf.rows();j++){
				for(int k=0;k<d0_buf.cols();k++)
					d0_buf(j,k) = d_temp0[j*d0_buf.cols()+k];
			}
		
		}
		
		//train step2
		
		if(party == ALICE){
			update_W0(W,tx_batch,txa_batch, b_IT2,d0, c_IT2, d1_buf); 
		}
		else{
			update_W1(W,tx_batch, txa_batch, b_IT2,d1, c_IT2, d0_buf);
		}
		
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
				cout<<res<<"\n";
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
	cout << "comm trash: " << (COMM_ / 512 / 1024) << "MB" << endl;
	cout << "\n";

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
