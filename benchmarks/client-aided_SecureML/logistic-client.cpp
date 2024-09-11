#include "utils.hpp"

double total_time, communication, computation, batchtime;

double offline;

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
    random_mat(a0);
    a1 = a - a0;
    send_mat(a0, io_client_alice);
    send_mat(a1, io_client_bob);
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
//    int t1 = mytime.start();

    Mat train_data(N,D), train_label(N,1);

    load_train_data(train_data, train_label);
    
//    cout<<mytime.end(t1)<<"s"<<endl;
    
    vector<int> perm = random_perm();
    
    io_client_alice->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_bob->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_alice->flush();
    io_client_bob->flush();

    sync_client();
    
    cout << "sending data......\n";

//    int t1 = mytime.start();

    Mat a(N, D);//, xa(N, D);
    random_mat(a);

    int t2 = mytime.start();

    client_random_distribute_mat(a);

//    xa = train_data - a;

    Mat b0(D, N * Ep), b1(B, IT);
    Mat c0(B, IT), c1(D, N * Ep);

    Mat b0_batch(D, B), b1_batch(B, 1);
    Mat c1_batch(D, B), c0_batch(B, 1);

    communication += mytime.end(t2);
    
    t2 = mytime.start();

    client_random_distribute_mat(train_data);
    client_random_distribute_mat(train_label);

    communication += mytime.end(t2);

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
    Mat a(N, D), xa(N, D);
    random_mat(a);
    xa = train_data - a;

    Mat b0(D, N * Ep), b1(B, IT);
    Mat c0(B, IT), c1(D, N * Ep);
    random_mat(b0);
    random_mat(b1);

    cout<<"calcing a,b,c......\n";

    for (int i = 0; i < IT; i++) {
        for (int j = 0; j < B; j++) {
            c0(j, i) = (a.row(perm[i * B + j]) * b0.col(i * B + j))(0, 0);
        //    cerr << "before c1\n";
            c1.col(i * B + j) = b1(j, i) * a.row(perm[i * B + j]).transpose();
        }
    }

    cout<<"sending a,b,c......\n";

    //client_random_distribute_mat(xa);
    send_mat(xa, io_client_alice);
    send_mat(xa, io_client_bob);
    client_random_distribute_mat(b0);
    client_random_distribute_mat(b1);
    client_random_distribute_mat(c0);
    client_random_distribute_mat(c1);*/
    
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
