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

#include <Eigen/Dense>

#define N 10240 //4096
#define D 785
#define B 128
#define testN 10000
#define Ep 2
#define IT N*Ep/B
#define L 20
#define P 64
#define EVE 3
#define EVALUATE 0
#define DEBUG_PRF 0
#define DEBUG_ACT 0

/*-------------------Malicious----------------*/
#define sigma 128
#define cN 32
#define ck 4
#define cM (ck * cN + ck)
/*-------------------Malicious----------------*/

/*-------------------Malicious-act----------------*/
#define dNN 8192
#define dk 3
/*-------------------Malicious-act----------------*/

/*-------------------Malicious-high-throughtput----------------*/
#define hN 1024
#define hk 5
#define	hM (hk * hN + hk + 1)
/*-------------------Malicious-high-throughtput----------------*/

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

double send_mat(Mat& x, NetIO * io) {
    vector<unsigned long int> temp(x.cols()*x.rows());
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            temp[j*x.cols()+k] = x(j,k);
        }
    }
    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
    io->flush();
    return sizeof(unsigned long int) * temp.size();
}

double receive_mat(Mat& x, NetIO * io) {
    vector<unsigned long int> temp(x.cols()*x.rows());
    io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
    for(int j=0; j<x.rows(); j++) {
        for(int k=0; k<x.cols(); k++) {
            x(j,k) = temp[j*x.cols()+k];
        }
    }
    return sizeof(unsigned long int) * temp.size();
}

void next_batch(Mat& batch,int start, vector<int>& perm, Mat& A) {
    for(int i=0; i<B; i++) {
        batch.row(i) = A.row(perm[start+i]);
    }
    return ;
}

Mat reconstruct(Mat A, int party, NetIO* io_server){
	Mat A_(A.rows(),A.cols());
	if(party==ALICE){
		send_mat(A, io_server);
		receive_mat(A_, io_server);
	}
	else{
		receive_mat(A_, io_server);
		send_mat(A, io_server);
	}
	return A + A_;
}

struct PRF{
	block key;
	unsigned long int index;
	
	PRF() {
		random_key();
		index = 0;
	}
	
	void random_key() {
		((unsigned long int*)&key)[0] = randomlong();
		((unsigned long int*)&key)[1] = randomlong();
		index = 0;
	}

	void cal(Mat& x) {
		AES_KEY k;
		AES_set_encrypt_key(key, &k);
		int r = x.rows(), c = x.cols();
		unsigned long int sz = ((unsigned long int)r * c + 1) >> 1;
		block *blks = new block[sz];
		for (unsigned long int i = 0; i < sz; ++i) {
			((unsigned long int*)&(blks[i]))[0] = 0;
			((unsigned long int*)&(blks[i]))[1] = i + index;
		}
		AES_ecb_encrypt_blks(blks, sz, &k);
		for (int i = 0, id = 0; i < r; ++i) for (int j = 0; j < c; ++j) {
			if((id & 1) == 0) x(i, j) = ((unsigned long int*)&(blks[id >> 1]))[0];
			else x(i, j) = ((unsigned long int*)&(blks[id >> 1]))[1];
			++id;
		}
		delete[] blks;
		index += sz;
	}
};

double send_keys(PRF *prf, NetIO *io, int sz) {
	vector<unsigned long int> temp(sz << 1);
	for (int i = 0; i < sz; ++i) {
		temp[i << 1] = ((unsigned long int*)&(prf[i].key))[0];
		temp[(i << 1) | 1] = ((unsigned long int*)&(prf[i].key))[1];
	}
	io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
	io->flush();
    return sizeof(unsigned long int) * temp.size();
}

double receive_keys(PRF *prf, NetIO *io, int sz) {
	vector<unsigned long int> temp(sz << 1);
	io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
	for (int i = 0; i < sz; ++i) {
		((unsigned long int*)&(prf[i].key))[0] = temp[i << 1];
		((unsigned long int*)&(prf[i].key))[1] = temp[(i << 1) | 1];
		prf[i].index = 0;
	}
    return sizeof(unsigned long int) * temp.size();
}


double send_half_mat(Mat& x, PRF *prf, NetIO *io, int o) {
	int r = x.rows(), c = x.cols();
    vector<unsigned long int> temp((r >> 1) * c);
    Mat y(1, c);
    for(int j = 0; j < (r >> 1); j++) {
    	prf[j + (r >> 1) * o].cal(y);
        for(int k = 0; k < c; k++) {
            temp[j * c + k] = x(j + (r >> 1) * o, k) - y(0, k);
        }
    }
    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
    io->flush();
    return sizeof(unsigned long int) * temp.size();
}

double receive_half_mat(Mat& x, PRF *prf, NetIO *io, int o) {
	int r = x.rows(), c = x.cols();
    vector<unsigned long int> temp((r >> 1) * c);
    io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
    Mat y(1, c);
    for(int j = 0; j < (r >> 1); j++) {
    	prf[j + (r >> 1) * o].cal(y);
    	x.row(j + (r >> 1) * o) = y;
        for(int k = 0; k < c; k++) {
            x(j + (r >> 1) * (1 - o), k) = temp[j * c + k];
        }
    }
    return sizeof(unsigned long int) * temp.size();
}

double distribute_mat(Mat& x, PRF *prf_alice, PRF *prf_bob, NetIO *io_alice, NetIO *io_bob) {
	assert((x.rows() & 1) == 0);
	double rs = 0;
	rs += send_half_mat(x, prf_alice, io_bob, 0);
	rs += send_half_mat(x, prf_bob, io_alice, 1);
	return rs;
}

double merge_mat(Mat& x, PRF *prf, NetIO *io_client, int party) {
	assert((x.rows() & 1) == 0);
	return receive_half_mat(x, prf, io_client, party - 1);
}

double send_half_mat_same(Mat& x, PRF *prf, NetIO *io, int o) {
	int r = x.rows(), c = x.cols();
    vector<unsigned long int> temp((r >> 1) * c);
    Mat y(1, c);
    for(int j = 0; j < (r >> 1); j++) {
    	prf->cal(y);
        for(int k = 0; k < c; k++) {
            temp[j * c + k] = x(j + (r >> 1) * o, k) - y(0, k);
        }
    }
    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
    io->flush();
    return sizeof(unsigned long int) * temp.size();
}

double receive_half_mat_same(Mat& x, PRF *prf, NetIO *io, int o) {
	int r = x.rows(), c = x.cols();
    vector<unsigned long int> temp((r >> 1) * c);
    io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
    Mat y(1, c);
    for(int j = 0; j < (r >> 1); j++) {
    	prf->cal(y);
    	x.row(j + (r >> 1) * o) = y;
        for(int k = 0; k < c; k++) {
            x(j + (r >> 1) * (1 - o), k) = temp[j * c + k];
        }
    }
    return sizeof(unsigned long int) * temp.size();
}

double distribute_mat_same(Mat& x, PRF *prf_alice, PRF *prf_bob, NetIO *io_alice, NetIO *io_bob) {
	assert((x.rows() & 1) == 0);
	double rs = 0;
	rs += send_half_mat_same(x, prf_alice, io_bob, 0);
	rs += send_half_mat_same(x, prf_bob, io_alice, 1);
	return rs;
}

double merge_mat_same(Mat& x, PRF *prf, NetIO *io_client, int party) {
	assert((x.rows() & 1) == 0);
	return receive_half_mat_same(x, prf, io_client, party - 1);
}

void add_mat(PRF *prf_alice, PRF *prf_bob, Mat& x) {
	Mat temp(x.rows(), x.cols());
	prf_alice->cal(x);
	prf_bob->cal(temp);
	x += temp;
}

void add_mat2(PRF *prf_alice, PRF *prf_bob, Mat& x) {
	assert(x.rows() == B);
	Mat temp(1, x.cols());
	for (int i = 0; i < B; ++i) {
		prf_alice[i].cal(temp);
		x.row(i) = temp;
		prf_bob[i].cal(temp);
		x.row(i) += temp;
	}
}

Mat sp_product(Mat x, Mat y) {
	Mat ret = x.row(0);
	for (int j = 0; j < x.cols(); j++) {
		ret(0, j) = 0;
		for (int i = 0; i < x.rows(); i++)
			ret(0, j) += x(i, j) * y(0, i);
	}
//	cerr << "sp_product: " << ret.rows() << ' ' << ret.cols() << endl;
	return ret;
}

#endif
