#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>

#include <time.h>

#include <emp-tool/emp-tool.h>
#include <emp-sh2pc/emp-sh2pc.h>

#include <Eigen/Dense>

#define N 10240
#define D 785
#define B 128
#define testN 10000
#define Ep 2
#define IT N*Ep/B
#define L 20
#define P 64
#define LAM 128
#define EVE 3
#define EVALUATE 0
#define DEBUG_PRF 0
#define DEBUG_ACT 0
#define DEBUG_TC 0

using namespace std;
using namespace Eigen;

typedef unsigned long int uint64;
typedef Matrix<unsigned long int,Dynamic,Dynamic> Mat;

struct Time_Calculator{
	#define CNT 10
	timeval tm[CNT];
	bool used[CNT];
	Time_Calculator() {
		for (int i = 0; i < CNT; ++i) used[i] = 0;
	}
	int start() {
		int o = -1;
		for (int i = 0; i < CNT; ++i) if(!used[i]) {
			used[i] = 1; o = i;
			break;
		}
		if(o == -1) {
			cout << "!!!!!!!!!!!!!!!!!!!!!\n";
			assert(0);
		}
		gettimeofday(&tm[o], NULL);
		return o;
	}
	double end(int o) {
		if(!used[o]) assert(0);
		timeval tmp;
		gettimeofday(&tmp, NULL);
		double rs = (tmp.tv_sec-tm[o].tv_sec) + (tmp.tv_usec-tm[o].tv_usec)/1000000.0;
		used[o] = 0;
		return rs;
	}
};

uint64 randomlong() {
    uint64 rand1 = abs(rand());
    uint64 rand2 = abs(rand());
    rand1 = rand1 << (sizeof(int)*8);
    uint64 randULL = (rand1 | rand2);
    return randULL;
}

void random_mat(Mat &x) {
    for (int i = 0; i < x.rows(); i++)
        for (int j = 0; j < x.cols(); j++)
            x(i, j) = randomlong();
}

int myrandom (int i) {
    return rand()%i;
}

vector<int> random_perm() {
    vector<int> temp,perm;
    for(int i=0; i<N; i++)
        temp.push_back(i);

    for(int i = 0; i<Ep; i++) {
        random_shuffle(temp.begin(),temp.end(),myrandom);
        perm.insert(perm.end(),temp.begin(),temp.end());

    }
    return perm;
}

void next_batch(Mat& batch,int start, vector<int>& perm, Mat& A) {
    for(int i=0; i<B; i++) {
        batch.row(i) = A.row(perm[start+i]);
    }
    return ;
}

void setup(int party, NetIO* &io_server, NetIO* &io_client_alice, NetIO* &io_client_bob, NetIO* &io_client) {
    if (party == ALICE || party == BOB) {
        io_server = new NetIO(party==ALICE ? nullptr : "127.0.0.1", 1025);
        io_server->set_nodelay();
    }

    if (party == ALICE || party == EVE) {
        io_client_alice = new NetIO(party==EVE ? nullptr : "127.0.0.1", 1026);
        io_client_alice->set_nodelay();
    }
    if (party == BOB || party == EVE) {
        io_client_bob = new NetIO(party==EVE ? nullptr : "127.0.0.1", 1027);
        io_client_bob->set_nodelay();
    }

    if(party == ALICE) io_client = io_client_alice;
    if(party == BOB) io_client = io_client_bob;
}

void send_mat(Mat& x, NetIO * io) {
    vector<unsigned long int> temp(x.cols()*x.rows());
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            temp[j*x.cols()+k] = x(j,k);
        }
    }
    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
    io->flush();
}

void send_mat(Mat& x, NetIO * io, double &COMM) {
    vector<unsigned long int> temp(x.cols()*x.rows());
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            temp[j*x.cols()+k] = x(j,k);
        }
    }
    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
    io->flush();

    COMM += (double)sizeof(unsigned long int) * temp.size();
}

void send_mat(Mat& x, NetIO *io, double& comm, double& comp, Time_Calculator& mytime) {
    
    int t1 = mytime.start();

    vector<unsigned long int> temp(x.cols()*x.rows());
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            temp[j*x.cols()+k] = x(j,k);
        }
    }

    comp += mytime.end(t1);
    
    t1 = mytime.start();

    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
    io->flush();

    comm += mytime.end(t1);
}

void receive_mat(Mat& x, NetIO * io) {
    vector<unsigned long int> temp(x.cols()*x.rows());
    io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            x(j,k) = temp[j*x.cols()+k];
        }
    }
}

void receive_mat(Mat& x, NetIO * io, double &COMM) {
    vector<unsigned long int> temp(x.cols()*x.rows());
    io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
    COMM += (double)sizeof(unsigned long int) * temp.size();
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            x(j,k) = temp[j*x.cols()+k];
        }
    }
}

#endif
