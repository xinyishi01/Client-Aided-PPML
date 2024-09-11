#include "utils.hpp"

double communication, computation, total_time;

Time_Calculator mytime;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

void load_train_data(Mat& train_data, Mat& train_label){
    ifstream infile( "../data/mnist_train.csv" );
        int count1=0, count2=0;
        int i=0;
        while(infile) {

            string s;
            if (!getline(infile,s))
                break;
            istringstream ss(s);
            int temp;
            char c;


            //read label
            ss>>temp;
            ss>>c;
            if(temp == 0 && count1<N/2) {
                train_label(i) = 0;
                count1++;

                //read data (last entry 1)
                for(int j=0; j<D-1; j++) {
                    ss>>train_data(i,j);
                    ss>>c;
                }

                train_data(i,D-1) = 1;
                i++;
            }


            if(temp != 0 && count2<N/2) {
                train_label(i) = 1;
                count2++;

                //read data (last entry 1)
                for(int j=0; j<D-1; j++) {
                    ss>>train_data(i,j);
                    ss>>c;
                }

                train_data(i,D-1) = 1;
                i++;
            }


            if(i>=N)
                break;
        }

        //train_data.conservativeResize(i, D);
        //train_label.conservativeResize(i,1);
    infile.close();
}

void client_random_distribute_mat(Mat a) {
    Mat a0, a1;
    a0 = a1 = a;

    int t1 = mytime.start();
    
    random_mat(a0);
    a1 = a - a0;

    computation += mytime.end(t1);

    send_mat(a0, io_client_alice, communication, computation, mytime);
    send_mat(a1, io_client_bob, communication, computation, mytime);
}

void sync_client() {
//    cerr << "start sync_client\n";
    int t1 = 666, t2;
    io_client_alice->recv_data(&t2, sizeof(int));
    io_client_bob->recv_data(&t2, sizeof(int));
    io_client_alice->send_data(&t1, sizeof(int));
    io_client_bob->send_data(&t1, sizeof(int));
    io_client_alice->flush();
    io_client_bob->flush();
//    cerr << "end sync_client\n";
}

int main(int argc, char** argv) {
    srand((unsigned)time(NULL));

    int port, party;
    parse_party_and_port(argv, &party, &port);
    setup(party, io_server, io_client_alice, io_client_bob, io_client);

    srand ( unsigned ( time(NULL) ) );
    
    total_time = 0.0;
    communication = 0.0;
    computation = 0.0;

    cout << "reading data......\n";

    Mat train_data(N,D), train_label(N,1);

    load_train_data(train_data, train_label);
    
    vector<int> perm = random_perm();
    
    io_client_alice->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_bob->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_alice->flush();
    io_client_bob->flush();
    
    sync_client();

    cout << "sending data......\n";

    int t1 = mytime.start();

    Mat a(N, D);//, xa(N, D);
    random_mat(a);
    client_random_distribute_mat(a);
//    xa = train_data - a;

    client_random_distribute_mat(train_data);
    client_random_distribute_mat(train_label);

    Mat b0(D, N * Ep), b1(B, IT);
    Mat c0(B, IT), c1(D, N * Ep);

    Mat b0_batch(D, B), b1_batch(B, 1);
    Mat c1_batch(D, B), c0_batch(B, 1);

    sync_client();

    for (int i = 0; i < IT; i++) {
        sync_client();
        random_mat(b0_batch);
        random_mat(b1_batch);
        for (int j = 0; j < B; j++) {
            c0_batch(j, 0) = (a.row(perm[i * B + j]) * b0_batch.col(j))(0, 0);
            c1_batch.col(j) = b1_batch(j, 0) * a.row(perm[i * B + j]).transpose();
        }
        client_random_distribute_mat(b0_batch);
        client_random_distribute_mat(b1_batch);
        client_random_distribute_mat(c0_batch);
        client_random_distribute_mat(c1_batch);
    }

/*
    random_mat(b0);
    random_mat(b1);

//    cout<<"calcing a,b,c......\n";

    for (int i = 0; i < IT; i++) {
        for (int j = 0; j < B; j++) {
            c0(j, i) = (a.row(perm[i * B + j]) * b0.col(i * B + j))(0, 0);
        //    cerr << "before c1\n";
            c1.col(i * B + j) = b1(j, i) * a.row(perm[i * B + j]).transpose();
        }
    }
*/
    computation += mytime.end(t1);

/*
    Mat b0(D, N), b1(1, N);
    Mat c0(1, N), c1(D, N);
    random_mat(b0);
    random_mat(b1);

    cout<<"calcing a,b,c......\n";

    for (int i = 0; i < N; i++) {
        c0(0, i) = (a.row(i) * b0.col(i))(0, 0);
        //    cerr << "before c1\n";
        c1.col(i) = b1(0, i) * a.row(i).transpose();
    }*/

    //cout<<"sending a,b,c......\n";

    //client_random_distribute_mat(xa);
    /*
    send_mat(xa, io_client_alice, communication, computation, mytime);
    send_mat(xa, io_client_bob, communication, computation, mytime);
    client_random_distribute_mat(b0);
    client_random_distribute_mat(b1);
    client_random_distribute_mat(c0);
    client_random_distribute_mat(c1);*/

    /*	
	PRF *prf_alice_temp = new PRF[N];
	PRF *prf_bob_temp = new PRF[N];
	
	receive_keys(prf_alice_temp, io_client_alice, N);
	receive_keys(prf_bob_temp, io_client_bob, N);
	
	distribute_mat(train_data, prf_alice_temp, prf_bob_temp, io_client_alice, io_client_bob);
	distribute_mat(train_label, prf_alice_temp, prf_bob_temp, io_client_alice, io_client_bob);
	
	delete [] prf_alice_temp;
	delete [] prf_bob_temp;

    int start = 0;
    Mat x_batch(B, D), y_batch(B, 1), tau0_batch(B, 1), tau1_batch(B, D);
    Mat z0(1, D), z1(B, 1);
    PRF *prf_alice = new PRF[B], *prf_alice_same = new PRF;
    PRF *prf_bob = new PRF[B], *prf_bob_same = new PRF;
    
    for (int i = 0; i < IT; ++i, start += B) {
    
    	next_batch(x_batch,start,perm,train_data);
        next_batch(y_batch,start,perm,train_label);
        
        //cout << "receive keys...\n";
        receive_keys(prf_alice, io_client_alice, B);
        receive_keys(prf_bob, io_client_bob, B);
        
        receive_keys(prf_alice_same, io_client_alice, 1);
        receive_keys(prf_bob_same, io_client_bob, 1);
        
        //cout << "add mat...\n";
        
        add_mat(prf_alice_same, prf_bob_same, z0);
        add_mat2(prf_alice, prf_bob, z1);
        
        
        //cout << "calculation...\n";
        
        for (int j = 0; j < B; ++j) {
        	tau0_batch(j, 0) = ((x_batch.row(j)).cwiseProduct(z0)).sum();
        	//tau1_batch.row(j) += x_batch.row(j) * z1(j, 0);
        	for (int k = 0; k < D; ++k) tau1_batch(j, k) = x_batch(j, k) * z1(j, 0);
        }
        
        //cout << "distribute mat...\n";
       
		distribute_mat(tau0_batch, prf_alice, prf_bob, io_client_alice, io_client_bob);
		distribute_mat(tau1_batch, prf_alice, prf_bob, io_client_alice, io_client_bob);
		
		if(DEBUG_PRF) {
			send_mat(tau1_batch, io_client_alice);
		}
        
    }

	delete[] prf_alice;
	delete[] prf_bob;
	delete prf_alice_same;
	delete prf_bob_same;*/
    
    
    //cout<<"total time:"<<total_time<<"s!!"<<endl;
    //cout<<"communication time:"<<communication<<"s!!"<<endl;
    //cout<<"computation time:"<<computation<<"s!!"<<endl;

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
