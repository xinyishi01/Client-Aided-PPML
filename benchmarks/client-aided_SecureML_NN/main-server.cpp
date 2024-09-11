//#pragma once
#include "defines.hpp"
#include "NeuralNetwork.hpp"

double total_time, online_time, offline_time, communication, computation, com_client;
double COMM, COMM_client_offline, COMM_client_online;

Time_Calculater mytime;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

SHOTIterated* ot;
SHOTIterated* ot2;

Mat train_data(N,D), train_label(N,10);
Mat x_batch(B,D), y_batch(B,10);
Mat test_data(testN,D), test_label(testN,1);
Mat aa(N, D), a[SIT][5][5], b[SIT][5][5], c[SIT][5][5], a_batch[5][5], b_batch[B][5][5], c_batch[B][5][5];
Mat xa(N, D), xa_batch(B, D), aa_batch(B, D);
unsigned long rs, r;

int start;

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

void load_test_data(Mat& test_data, Mat& test_label)
{
    ifstream infile2( "../data/mnist_test.csv" );

    int i=0;

    cout<<"load testing data.......\n";

    while(infile2) {

        string s;
        if (!getline(infile2,s)) {
            break;
        }
        istringstream ss(s);
        int temp;
        char c;

        //read label
        ss>>temp;
        ss>>c;

        //if(temp == 0 || temp == 1){
        test_label(i) = temp;


        //read data (last entry 1)
        for(int j=0; j<D; j++) {
            ss>>test_data(i,j);
            test_data(i,j) <<= L;
            ss>>c;
        }

    //    test_data(i,D-1) = 1;
        i++;
        //}

    }

    test_data.conservativeResize(i, D);
    test_label.conservativeResize(i,1);

    infile2.close();

    return;
}

void sync_server() {
//  cerr << "start sync_server\n";
    int t1 = 666, t2;
    io_client->send_data(&t1, sizeof(int));
    io_client->flush();
    io_client->recv_data(&t2, sizeof(int));
//  cerr << "end sync_server\n";
}

void sync_between_servers(int party) {
//    cerr << "start sync_bet_server\n";
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
//    cerr << "end sync_bet_server\n";
}

vector<int> prework(NeuralNetwork* &net, int party, int port) {
    //setup connection
    setup_semi_honest(io_server, party);

    ot = new SHOTIterated(io_server, party == ALICE);
    ot2 = new SHOTIterated(io_server, party != ALICE);

    net = new NeuralNetwork(party, io_server);

    for (int i = 1; i <= LL; ++i) {
    //    cout << i << endl;
        net->layers[i]->set_io(io_server, ot, ot2);
    }

    vector<int> perm(N * Ep, 0);

    io_client->recv_data(&perm[0], sizeof(int) * perm.size());

    sync_between_servers(party);
    sync_server();

    total_time = 0.0;
    communication = 0.0;
    computation = 0.0;
    com_client = 0.0;

    COMM = 0.0;
    COMM_client_online = 0.0;
    COMM_client_offline = 0.0;

    cout << "before load_train_data\n";
    
    int t1 = mytime.start();

    receive_mat(train_data, io_client, COMM_client_online);
    receive_mat(train_label, io_client, COMM_client_online);

    for (int i = 0; i < train_data.rows(); i++)
        for (int j = 0; j < train_data.cols(); j++)
            train_data(i, j) <<= L;
    for (int i = 0; i < train_label.rows(); i++)
        for (int j = 0; j < train_label.cols(); j++)
            train_label(i, j) <<= L;

    online_time += mytime.end(t1);

 //   for (int i = 0; i < train_data.rows(); i++) {
    /*
        for (int j = 0; j < train_data.cols(); j++)
            cerr << reconstruct(train_data, party)(0, j) << ' ';
        cerr << '\n';*/

    if(party == ALICE) {
        load_test_data(test_data, test_label);
    }
/*
    cout << "before read_MT\n";

    vector<int> perm = read_MT(party-1);

    cout << "after read_MT\n";*/

 //   xa = reconstruct(train_data - aa, party);

    net->layers[1]->inputDim = D;
    net->layers[1]->outputDim = net->layers[2]->inputDim = LAYER1;
    net->layers[2]->outputDim = net->layers[3]->inputDim = LAYER2;
    net->layers[3]->outputDim = LAST_LAYER_SIZE;

    cout << "before initialize\n";

    net->layers[1]->initialize(LearningRate);
    net->layers[2]->initialize(LearningRate);
    net->layers[3]->initialize(LearningRate);
/*
    for (int i = 1; i <= LL; ++i) {
    //    cout << i << endl;
        net->layers[i]->initialize();
    }*/

    cout << "after initialize\n";

    int lay[4];
    lay[0] = D;
    lay[1] = LAYER1;
    lay[2] = LAYER2;
    lay[3] = 10;
    for (int j = 1; j <= LL; j++) {
        a_batch[j][1].resize(B, lay[j - 1]);
        for (int k = 0; k < B; k++)
            b_batch[k][j][1].resize(lay[j - 1], lay[j]);
        c_batch[0][j][1].resize(B, lay[j]);
                
        a_batch[j][2].resize(B, lay[j]);
        b_batch[0][j][2].resize(B, lay[j]);
        c_batch[0][j][2].resize(B, lay[j]);

        a_batch[j][3].resize(B, lay[j]);
        for (int k = 0; k < B; k++)
            b_batch[k][j][3].resize(lay[j], lay[j - 1]);
        c_batch[0][j][3].resize(B, lay[j - 1]);

        a_batch[j][4].resize(lay[j - 1], B);
        b_batch[0][j][4].resize(B, lay[j]);
        for (int k = 0; k < B; k++)
            c_batch[k][j][4].resize(lay[j - 1], lay[j]);
    }

    return perm;
}

void readMiniBatch(NeuralNetwork* net, int d, string phase, vector<int>& perm, int party) {
 //   cout << "begin\n";
    next_batch(x_batch, start, perm, train_data);
    next_batch(y_batch, start, perm, train_label);
    net->inputData = x_batch;
    net->outputData = y_batch;/*
    next_batch(xa_batch, B * d, perm, xa);
    next_batch(aa_batch, B * d, perm, aa);*/

    int t1 = mytime.start();
    int t2 = mytime.start();

    for (int i = 1; i <= LL; ++i)
        for (int j = 1; j < 5; ++j) {
            /*
            a_batch[i][j] = a[d][i][j];
            b_batch[i][j] = b[d][i][j];
            c_batch[i][j] = c[d][i][j];*/
            receive_mat(a_batch[i][j], io_client, COMM_client_offline);
            if (j == 2 || j == 4)
                receive_mat(b_batch[0][i][j], io_client, COMM_client_offline);
            else {
                for (int o = 0; o < B; o++)
                    receive_mat(b_batch[o][i][j], io_client, COMM_client_offline);
            }
            if (j != 4)
                receive_mat(c_batch[0][i][j], io_client, COMM_client_offline);
            else {
                for (int o = 0; o < B; o++)
                    receive_mat(c_batch[o][i][j], io_client, COMM_client_offline);
            }
        }
    io_client->recv_data(&rs, sizeof(unsigned long));
    io_client->recv_data(&r, sizeof(unsigned long));

    offline_time += mytime.end(t1);
    total_time += mytime.end(t2);

    COMM_client_offline += 2 * sizeof(unsigned long);

    int p = 1;
    Integer rs_ = Integer(64, rs, p);
    net->layers[LL]->set_r(rs_, r);
//    cerr << rs << ' ' << r << endl;

    net->layers[1]->is_first = 0;
//    net->layers[1]->set_xa(xa_batch);
//    cerr << "before set_MT\n";
    for (int i = 1; i <= LL; ++i) {
    //    cerr << a_batch[i][1].rows() << ' ' << a_batch[i][1].cols() << endl;
        net->layers[i]->a1 = a_batch[i][1];
        net->layers[i]->a2 = a_batch[i][2];
        net->layers[i]->a3 = a_batch[i][3];
        net->layers[i]->a4 = a_batch[i][4];

        for (int o = 0; o < B; o++)
            net->layers[i]->b1[o] = b_batch[o][i][1];
        net->layers[i]->b2 = b_batch[0][i][2];
        for (int o = 0; o < B; o++)
            net->layers[i]->b3[o] = b_batch[o][i][3];
        net->layers[i]->b4 = b_batch[0][i][4];
        
        net->layers[i]->c1 = c_batch[0][i][1];
        net->layers[i]->c2 = c_batch[0][i][2];
        net->layers[i]->c3 = c_batch[0][i][3];
        net->layers[i]->c4 = c_batch[0][i][4];
        for (int o = 1; o < B; o++)
            net->layers[i]->c4 += c_batch[o][i][4];

    //    net->layers[i]->set_MT(a_batch[i][1], b_batch[i][1], c_batch[i][1], a_batch[i][2], b_batch[i][2], c_batch[i][2], a_batch[i][3], b_batch[i][3], c_batch[i][3], a_batch[i][4], b_batch[i][4], c_batch[i][4]);
    }
//    cerr << "after set_MT\n";
    start += B;
}

void truncation(myType& x) {
    x = (long int)x >> L;
}

void trunc(Mat& x) {
    for (int i = 0; i < x.rows(); ++i)
        for (int j = 0; j < x.cols(); ++j) {
            truncation(x(i, j));
        }
}

Mat plain_ReLU(Mat& x, Mat& W, Mat& bias) {
    Mat rs = x * W;
    for (int i = 0; i < rs.rows(); ++i) rs.row(i) += bias;
    for (int i = 0; i < rs.rows(); ++i) for (int j = 0; j < rs.cols(); ++j) {
        if((long int)rs(i, j) < 0) rs(i, j) = 0;
    }
    trunc(rs);
    return rs;
}

void train(NeuralNetwork* net, vector<int>& perm, int party) {
//    cout << "train\n";

    start = 0;

    double mx = 0;

    sync_between_servers(party);
    sync_server();

    total_time = online_time + offline_time;

    for (int i = 0; i < IT; ++i) {
        // cout << "----------------------------------" << endl;  
        cout << "Iteration " << i << endl;

        sync_server();
        
        readMiniBatch(net, i % SIT, "TRAINING", perm, party);

    //    cout << "after mini\n";

        int t1 = mytime.start();
        int t2 = mytime.start();

     //   int t3 = mytime.start();

        net->forward();
    //    cerr << "after forward\n";
    //    cerr << "total forward" 

        // start_m();
        net->backward();

        // cout << "----------------------------------" << endl;

        online_time += mytime.end(t1);
        total_time += mytime.end(t2);
        
        if(EVALUATE && i % 5 == 0) {
            //cout << "evaluate ...\n";
            Mat W1 = net->layers[1]->reconstruct(net->layers[1]->W);
            Mat bias1 = net->layers[1]->reconstruct(net->layers[1]->bias);
            Mat W2 = net->layers[2]->reconstruct(net->layers[2]->W);
            Mat bias2 = net->layers[2]->reconstruct(net->layers[2]->bias);
            Mat W3 = net->layers[3]->reconstruct(net->layers[3]->W);
            Mat bias3 = net->layers[3]->reconstruct(net->layers[3]->bias);
            if(party == ALICE) {
            /*    Mat X1 = plain_ReLU(test_data, W1, bias1);
                Mat X2 = plain_ReLU(X1, W2, bias2);
                Mat X3 = plain_ReLU(X2, W3, bias3);
                */
                Mat X1 = plain_ReLU(test_data, W1, bias1);
                Mat X2 = plain_ReLU(X1, W2, bias2);
                Mat X3 = plain_ReLU(X2, W3, bias3);

                double tt = 0, hs = 0.001;
                for(int i = 0; i < X3.rows(); ++i) for (int j = 0; j < X3.cols(); ++j) {
                    if((long int)X3(i, j) > 0) {
                        tt += X3(i, j);
                        ++hs;
                    }
                }
                tt /= pow(2, L);


                int ac = 0;
                for (int i = 0; i < testN; ++i) {
                    unsigned long int maxvalue = X3(i, 0);
                    int label = 0;
                    for (int j = 1; j < 10; ++j) {
                        if(X3(i, j) > maxvalue) {
                            maxvalue = X3(i, j);
                            label = j;
                        }
                    }
                    if(label == test_label(i, 0)) ++ac;
                //    cout << test_label(i, 0) << '&' << label << ' ';
                }
            //    cout << '\n';
                mx = max(mx, (double) 100.0 * ac / testN);
                cout << "#" << i <<" accuracy : " << (double) 100.0 * ac / testN << "%\n";
                cout << "#" << i <<" max accuracy : " << mx << "%\n";
                cerr << "................................" << tt / hs << endl;

            }
        }
    }

}

int main(int argc, char** argv) {
    int port, party;
    parse_party_and_port(argv, &party, &port);
    setup(party, io_server, io_client_alice, io_client_bob, io_client);
    NeuralNetwork* network;
    vector<int> perm = prework(network, party, port);
    cout << "prework finished\n";
    train(network, perm, party);

    cout << "total time:" << total_time << "s" << endl;
    cout << "online time:" << online_time << "s" << endl;
    cout << "offline_time:" << offline_time << "s" << endl;
    cout << "comm between servers: " << (COMM / 512 / 1024) << "MB" << endl;
    cout << "online comm between servers and clients: " << (COMM_client_online / 512 / 1024) << "MB" << endl;
    cout << "offline comm between servers and clients: " << (COMM_client_offline / 512 / 1024) << "MB" << endl;
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
