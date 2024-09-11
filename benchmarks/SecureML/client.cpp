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
    

    cout << "sending data......\n";
    
    Mat a(N, D);//, xa(N, D);
    random_mat(a);
    client_random_distribute_mat(a);
    
    sync_client();

    client_random_distribute_mat(train_data);
    client_random_distribute_mat(train_label);
	
	Mat a_batch(B, D);
    Mat b0_batch(D, 1), b1_batch(B, 1);
    Mat c0_batch, c1_batch;

    for (int i = 0; i < IT; i++) {
        random_mat(b0_batch);
        random_mat(b1_batch);
        for (int j = 0; j < B; ++j) {
        	a_batch.row(j) = a.row(perm[i * B + j]);
        }
        c0_batch = a_batch * b0_batch;
        c1_batch = a_batch.transpose() * b1_batch;
        client_random_distribute_mat(b0_batch);
        client_random_distribute_mat(b1_batch);
        client_random_distribute_mat(c0_batch);
        client_random_distribute_mat(c1_batch);
    }

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
