
#ifndef UTILS_HPP
#define UTILS_HPP


#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <set>

#include <time.h>

#include <emp-tool/emp-tool.h>
#include <emp-sh2pc/emp-sh2pc.h>

#include <Eigen/Dense>


#define N 10240
#define D 784
#define K 10
#define B 128
#define testN 10000
#define Ep 2
#define IT ((N*Ep)/B)
#define L 15
#define P 64
#define EVALUATE 0
#define DEBUG_FORWARD 0
#define DEBUG_BACKWARD 0
#define DEBUG_INNER 0
#define DEBUG_ACT 0
#define DEBUG_TC 0
#define learning_rate (1 << 7)

#define EVE 3


/*--------------------------Malicious--------------------------*/
#define cN 16//128
#define ck 4
#define cM (cN * ck + ck)

#define rN 16//128
#define rk 5
#define rM (rN * rk + rk + 1)

#define dNN 20480
#define dk 2
/*--------------------------Malicious--------------------------*/
	


using namespace std;
using namespace Eigen;


typedef unsigned long int uint64;
typedef Eigen::Matrix<unsigned long int,Eigen::Dynamic,Eigen::Dynamic> Mat;

typedef uint8_t* vec_t;

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
		if(!used[o]) {
			cout << "!!!!!!!!!!!!!!!!!!!!!\n";
			assert(0);
		}
		timeval tmp;
		gettimeofday(&tmp, NULL);
		double rs = (tmp.tv_sec-tm[o].tv_sec) + (tmp.tv_usec-tm[o].tv_usec)/1000000.0;
		used[o] = 0;
		return rs;
	}
};


uint64 randomlong(){
	uint64 rand1 = abs(rand());
    uint64 rand2 = abs(rand());
    rand1 = rand1 << (sizeof(int)*8);   
    uint64 randULL = (rand1 | rand2);   
    return randULL;
}

void rand_mat(Mat& x) {
	for (int i = 0; i < x.rows(); ++i) for (int j = 0; j < x.cols(); ++j) 
		x(i, j) = randomlong();
}

void debug_mat(Mat& x, string str) {
	ofstream ofs(str);
	for (int i = 0; i < x.rows(); ++i) {
 		for (int j = 0; j < x.cols(); ++j) {
 				ofs << (double)((long int)x(i, j)) / (double)pow(2, L);
 				ofs << " ";
  		}
  		ofs << "\n";
	}
	ofs.close();
}

uint64 myrand() {
	const unsigned long O = 100, H = 100;
	unsigned long int rd = (rand() * rand()) % (O + H);
	rd = rd * (1 << L) / 10000;
	return rd - H * (1 << L) / 10000;
}

int myrandom (int i) { return rand()%i;}

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

void next_batch(Mat& batch,int start, vector<int>& perm, Mat& A) {
    for(int i = 0; i < B; ++i) {
        batch.row(i) = A.row(perm[start+i]);
    }
    return ;
}

//#define trunc(x) (x = ((long int)x) >> L)

void trunc(uint64& x) {
	x = (long int)x >> L;
}

void truncation_mat(Mat& x) {
	for (int i = 0; i < x.rows(); ++i) for (int j = 0; j < x.cols(); ++j) {
		trunc(x(i, j));
	}
}

Mat plain_ReLU(const Mat& x, Mat& W, Mat& bias) {
	Mat rs = x * W;
	truncation_mat(rs);
	for (int i = 0; i < rs.rows(); ++i) rs.row(i) += bias;
	for (int i = 0; i < rs.rows(); ++i) for (int j = 0; j < rs.cols(); ++j) {
		if((long int)rs(i, j) < 0) rs(i, j) = 0;
	}
	return rs;
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

Mat flatten_mat(const Mat& x) {
	int r = x.rows(), c = x.cols();
	Mat rs(1, r * c);
	for (int i = 0, id = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j, ++id) rs(id) = x(i, j);
	}
	return rs;
}

void unfold_mat(Mat& x, Mat y) {
	int r = x.rows(), c = x.cols();
	for (int i = 0, id = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j, ++id) x(i, j) = y(id);
	}
}

void add_mat(PRF *prf_alice, PRF *prf_bob, Mat& x) {
	Mat temp(x.rows(), x.cols());
	prf_alice->cal(x);
	prf_bob->cal(temp);
	x += temp;
}

void cal_mat2(PRF *prf, Mat& x) {
	assert(x.rows() == B);
	Mat temp(1, x.cols());
	for (int i = 0; i < B; ++i) {
		prf[i].cal(temp);
		x.row(i) = temp;
	}
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
#endif
