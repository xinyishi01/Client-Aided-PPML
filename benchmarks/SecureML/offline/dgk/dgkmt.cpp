#include <iostream>
//#include "pailler.h"
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include "dgk.h"

#include <time.h>

#include <Eigen/Dense>
#include <emp-tool/emp-tool.h>

#define N 1024 
//#define N 12665
#define D 785
#define B 128
#define Ep 2
#define IT (N*Ep/B)
#define P 64

using namespace std;
using namespace Eigen;

typedef unsigned long int uint64;
typedef Matrix<unsigned long int,Dynamic,Dynamic> Mat;

NetIO * io;

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
}mytime;

Mat reconstruct(Mat A, int party){
	vector<unsigned long int> A_temp(A.cols()*A.rows());
	
	int tempi = 10;
	if(party==ALICE){
		io->send_data(&tempi,4);
		io->recv_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
		Mat A_(A.rows(),A.cols());
		for(int i=0;i<A_.rows();i++){
			for(int j=0;j<A_.cols();j++){
				A_(i,j) = A_temp[i*A_.cols()+j];
			}
		}
		
		Mat A_rec = A+A_;
		
		/*
		for(int i=0;i<A_rec.rows();i++){
			for(int j=0;j<A_rec.cols();j++)
				cout<<A_rec(i,j)<<",";
			cout<<";";
		}
		cout<<endl;
		*/
		
		//cout<<A_rec<<endl;
		
		return A_rec;
	}
	else{
		for(int i=0;i<A.rows();i++){
			for(int j=0;j<A.cols();j++)
				A_temp[i*A.cols()+j] = A(i,j);
		}
		
		io->recv_data(&tempi,4);
		io->send_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
		
		
		
		return A;
	}


}


int myrandom (int i) { return rand()%i;}

uint64 randomlong(){
	uint64 rand1 = abs(rand());
    uint64 rand2 = abs(rand());
    rand1 = rand1 << (sizeof(int)*8);   
    uint64 randULL = (rand1 | rand2);   
    return randULL;

}

vector<int> random_perm(){
	vector<int> temp,perm;
	for(int i=0;i<N;i++)
		temp.push_back(i);
	
	for(int i = 0;i<Ep;i++){
		random_shuffle(temp.begin(),temp.end(),myrandom);
		perm.insert(perm.end(),temp.begin(),temp.end());
		
	}
	return perm;
}

void next_batch(Mat& batch,int start, vector<int>& perm, Mat& A){
	
	
	for(int i=0;i<B;i++){
		batch.row(i) = A.row(perm[start+i]);
	}
	return ;
}

void same_base_mul( mpz_t* c,  mpz_t b, unsigned long* a, int n,  int bitl, mpz_t p){
	mpz_t powers[bitl];
	for(int i=0;i<bitl;i++){
		mpz_init(powers[i]);
	}
	mpz_set(powers[0],b);
	for(int i=1;i<bitl;i++){
		//mpz_mul(powers[i],powers[i-1],powers[i-1]);
		//mpz_mod(powers[i],powers[i],p);
		mpz_powm_ui(powers[i],powers[i-1],2,p);
	}
	
	for(int i=0;i<n;i++){
		unsigned long temp = a[i];
		mpz_set_ui(c[i],1);
		for(int j=0;j<bitl;j++){
			if(temp&1){
				mpz_mul(c[i],c[i],powers[j]);
				mpz_mod(c[i],c[i],p);
			}
			temp>>=1;
		}
	}

	return;

}

int main(int argc, char** argv) {
	
	
	
	srand ( unsigned ( time(NULL) ) );
	
	PRG prg;
	
	gmp_randstate_t rnd;
	gmp_randinit_default(rnd);
	gmp_randseed_ui(rnd, rand());
	dgk_pubkey_t * pub;
	dgk_prvkey_t * prv;
	
	
	
	unsigned int l = 64;
	unsigned int nbit = 2048;
	
	dgk_keygen(nbit, l, &pub, &prv);
	//dgk_readkey(nbit, l, &pub, &prv);
	
	int port, party;
	parse_party_and_port(argv, &party, &port);
	io = new NetIO(party==ALICE ? nullptr : "127.0.0.1", port);
	io->set_nodelay();
	
	int t2 = mytime.start();
	Mat a(N,D),b1(D,IT),c1(B,IT),b2(B,IT),c2(D,IT);
	
	vector<int> perm = random_perm();
	
	
	
	for(int i=0;i<a.rows();i++){
		for(int j=0;j<a.cols();j++){
			a(i,j) = randomlong();
		}
	}
	
	for(int i=0;i<b1.rows();i++){
		for(int j=0;j<b1.cols();j++){
			b1(i,j) = randomlong();
		}
	}
	
	for(int i=0;i<b2.rows();i++){
		for(int j=0;j<b2.cols();j++){
			b2(i,j) = randomlong();
		}
	}
	
	
	
	
	int start_setup = 0;
	for(int round=0;round<IT;round++){
		Mat a_batch(B,D), a_batch_t(D,B);
		
		next_batch(a_batch, start_setup,perm,a);

		int t1=mytime.start();
		
		cout<<"enc:\n";
		
		mpz_t b_m[D],cipher[D];
		for(int i=0;i<D;i++){
			mpz_init(b_m[i]);
			mpz_init(cipher[i]);
			mpz_set_ui(b_m[i],b1(i, round));
			
			dgk_encrypt_crt(cipher[i], pub, prv, b_m[i], rnd);

			

		}
		
		cout<<mytime.end(t1)<<"s enc"<<endl;
		
		t1 = mytime.start();
		
		
		
		cout<<"exchange:\n";
		
		
		char temp[410*D],temp2[410*D];
		for(int i=0;i<D;i++){
			mpz_get_str((temp+i*410),62,cipher[i]);
		
		}
		
		
		
		if(party == ALICE){
			
			
			io->send_data(temp, 410*D);
			io->recv_data(temp2, 410*D);
		
		}
		else{
			io->recv_data(temp2, 410*D);
			io->send_data(temp, 410*D);
			
		}
		
		io->flush();
		
		for(int i=0;i<D;i++){
			mpz_set_str(cipher[i],(temp2+i*410),62);
		
		}
		
		/*
		if(party!=ALICE){
			for(int i=0;i<D;i++){
				
				cout<<mpz_out_str(stdout,62,cipher[i])<<","<<endl;
			
			}
		
		}
		*/
		
		
		cout<<mytime.end(t1)<<"s exchange"<<endl;
		
		t1 = mytime.start();
		cout<<"homo:\n";
		
		
		mpz_t cipher2[D][B];
		
		
		
		
		for(int i=0;i<D;i++){
			for(int j=0;j<B;j++){
				mpz_init(cipher2[i][j]);
			}
		}
		
		for(int i=0;i<D;i++){
			
			same_base_mul(cipher2[i],cipher[i],a_batch.col(i).data(),B,64,pub->n);
		}
		
		cout<<mytime.end(t1)<<"s mul"<<endl;
		t1=mytime.start();
		
		
		
		mpz_t cipher3[B];
		mpz_t r1[B];
		mpz_t r1_enc[B];
		
		for(int i=0;i<B;i++){
			mpz_init(cipher3[i]);
			mpz_set_ui(cipher3[i],1);
			
			mpz_init(r1[i]);
			//mpz_set_ui(r1[i],0);
			prg.random_mpz(r1[i], 64);
			//mpz_mod(r1[i],r1[i],pk2.n);
			
			mpz_init(r1_enc[i]);
			dgk_encrypt_crt(r1_enc[i], pub, prv, r1[i], rnd);
			//paillier.Enc(r1_enc[i],r1[i],pk2);
			
			for(int j=0;j<D;j++){
				mpz_mul(cipher3[i],cipher3[i],cipher2[j][i]);
				mpz_mod(cipher3[i],cipher3[i],pub->n);
			}
			
			mpz_mul(cipher3[i],cipher3[i],r1_enc[i]);
			mpz_mod(cipher3[i],cipher3[i],pub->n);
			
			//paillier.Add(cipher3[i],cipher3[i],r1_enc[i],pk2);
		}
		
		cout<<mytime.end(t1)<<"s add"<<endl;
		t1 = mytime.start();
		
		
		
		cout<<"exchange:\n";
		
		
		char temp3[410*B],temp4[410*B];
		for(int i=0;i<B;i++){
			mpz_get_str((temp3+i*410),62,cipher3[i]);
		
		}
		
		if(party == ALICE){
			io->send_data(temp3, 410*B);
			io->recv_data(temp4, 410*B);
		
		}
		else{
			io->recv_data(temp4, 410*B);
			io->send_data(temp3, 410*B);
		}
		
		io->flush();
		
		for(int i=0;i<B;i++){
			mpz_set_str(cipher3[i],(temp4+i*410),62);
		
		}
		
		cout<<mytime.end(t1)<<"s exchange"<<endl;
		t1 = mytime.start();
		
		cout<<"Dec:\n";
		
		mpz_t c_m[B];
		for(int i=0;i<B;i++){
			mpz_init(c_m[i]);
			dgk_decrypt(c_m[i], pub, prv, cipher3[i]);
			
			c1(i,round) = mpz_get_ui(c_m[i])-mpz_get_ui(r1[i]);
		}
		
		c1.col(round)+= a_batch*b1.col(round);
		
		cout<<mytime.end(t1)<<"s dec"<<endl;
		
		/*
		Mat a_vec,b_vec,c_vec;
		
		a_vec = reconstruct(a_batch,party);
		b_vec = reconstruct(b1.col(round),party);
		c_vec = reconstruct(c1.col(round),party);
		
		cout<<(a_vec*b_vec).isApprox(c_vec)<<endl;
		*/
		
		start_setup+=B;
		
	}
	
	cout<<"total time: "<<mytime.end(t2) << "s!!\n";
	
}
