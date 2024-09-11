#include "defines.hpp"

double communication, computation, total_time;

Time_Calculater mytime;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

Mat a0(N,D),a1(N,D);
Mat a[5][5], b[B][5][5], c[B][5][5];
Mat test_data(testN,D), test_label(testN,1);

void client_random_distribute_mat(Mat a) {
    Mat a0, a1;
    a0 = a1 = a;
    random_mat(a0);
    a1 = a - a0;
/*
    send_mat(a0, io_client_alice, communication, computation, mytime);
    send_mat(a1, io_client_bob, communication, computation, mytime);*/
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


void load_train_data(Mat& train_data, Mat& train_label) {
    ifstream infile( "../data/mnist_train.csv" );
    int i=0;
    cout<<"load train data.......\n";
    
    while(infile) {
    //  cout << i << endl;

        string s;
        if (!getline(infile,s))
            break;
        istringstream ss(s);
        int temp;
        char c;
        
        //read label
        ss>>temp;
        ss>>c;
        for (int j = 0; j < 10; j++)
            train_label(i, j) = 0;
        train_label(i, temp) = 1;//0~9
        
        //read data 
        for(int j=0; j<D; j++) {
            ss>>train_data(i,j);
            ss>>c;
        }
        i++;
        if(i>=N) break;
    }

    train_data.conservativeResize(i, D);
    train_label.conservativeResize(i, 10);

    cout<<"n = "<<i<<endl;
    infile.close();
}

int main(int argc, char** argv){
    srand((unsigned)time(NULL));

    int port, party;
    /*
    parse_party_and_port(argv, &party, &port);
    setup(party, io_server, io_client_alice, io_client_bob, io_client);*/

    srand(unsigned(time(NULL)));

    Mat train_data(N,D), train_label(N,10);
    
    load_train_data(train_data, train_label);

    vector<int> perm = random_perm();
    /*
    io_client_alice->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_bob->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_alice->flush();
    io_client_bob->flush();*/
    
//    sync_client();

    cout << "sending data......\n";
    client_random_distribute_mat(train_data);
    client_random_distribute_mat(train_label);

    total_time = 0;
    int t1 = mytime.start();

    int lay[4];
    lay[0] = D;
    lay[1] = LAYER1;
    lay[2] = LAYER2;
    lay[3] = 10;

//    sync_client();

    for (int i = 0; i < IT; i++) {
//        sync_client();
        cerr << i << endl;

        for (int j = 1; j <= LL; j++) {
        //    cerr << "!" << j << endl;
            a[j][1].resize(B, lay[j - 1]);
            for (int k = 0; k < B; k++)
                b[k][j][1].resize(lay[j - 1], lay[j]);
            c[0][j][1].resize(B, lay[j]);
                
            a[j][2].resize(B, lay[j]);
            b[0][j][2].resize(B, lay[j]);
            c[0][j][2].resize(B, lay[j]);

            a[j][3].resize(B, lay[j]);
            for (int k = 0; k < B; k++)
                b[k][j][3].resize(lay[j], lay[j - 1]);
            c[0][j][3].resize(B, lay[j - 1]);

            a[j][4].resize(lay[j - 1], B);
            b[0][j][4].resize(B, lay[j]);
            for (int k = 0; k < B; k++)
                c[k][j][4].resize(lay[j - 1], lay[j]);
        }
        for (int j = 1; j <= LL; j++)
            for (int k = 1; k <= 4; k++) {
                random_mat(a[j][k]);
                if (k == 2 || k == 4)
                    random_mat(b[0][j][k]);
                else {
                    for (int o = 0; o < B; o++)
                        random_mat(b[o][j][k]);
                }
            }

    //    cerr << "1\n";

        for (int j = 1; j <= LL; j++)
            for (int k = 1; k <= 4; k++) {
            //    cerr << k << endl;
                if (k == 2) {
                    for (int mi = 0; mi < c[0][j][k].rows(); mi++)
                        for (int mj = 0; mj < c[0][j][k].cols(); mj++)
                            c[0][j][k](mi, mj) = a[j][k](mi, mj) * b[0][j][k](mi, mj);
                } else
                if (k == 1 || k == 3) {
                    for (int o = 0; o < B; o++)
                        c[0][j][k].row(o) = (a[j][k].row(o)) * b[o][j][k];
                } else
                if (k == 4) {
                    for (int o = 0; o < B; o++) {
                        c[o][j][k] = a[j][k].col(o) * b[0][j][k].row(o);
                    }
                }
/*
                if (k != 2) {
                    c[j][k] = a[j][k] * b[j][k];
                } else
                if (k == 1) {

                }
                else {
                    
                }*/
            }

    //    cerr << "2\n";

        for (int j = 1; j <= LL; j++)
            for (int k = 1; k <= 4; k++) {
                client_random_distribute_mat(a[j][k]);
            //    client_random_distribute_mat(b[j][k]);
                if (k == 2 || k == 4)
                    client_random_distribute_mat(b[0][j][k]);
                else {
                    for (int o = 0; o < B; o++)
                        client_random_distribute_mat(b[o][j][k]);
                }
                if (k != 4)
                    client_random_distribute_mat(c[0][j][k]);
                else {
                    for (int o = 0; o < B; o++)
                        client_random_distribute_mat(c[o][j][k]);
                }
            //    client_random_distribute_mat(c[j][k]);
            }

        //------------division----------------
        /*
        bool t;
        unsigned long rs, r1, r2;
        rs = rand();
        r1 = r2 = 0;
        for (int j = 63; j >= 0; --j) {
            if ((rs >> j) & 1) {
                t = rand() & 1;
                if (t)
                    r2 += (1 << j);
                else
                    r1 += (1 << j);
            }
        }*/
        /*
        io_client_alice->send_data(&rs, sizeof(unsigned long));
        io_client_alice->send_data(&r1, sizeof(unsigned long));
        io_client_alice->flush();
        io_client_bob->send_data(&rs, sizeof(unsigned long));
        io_client_bob->send_data(&r2, sizeof(unsigned long));
        io_client_bob->flush();*/
    }

    total_time = mytime.end(t1);
    cerr << "total_time:" << total_time << "s\n";
/*
    if (party != EVE) {
        delete io_server;
    }
    if (party != BOB) {
        delete io_client_alice;
    }
    if (party != ALICE) {
        delete io_client_bob;
    }*/

}
