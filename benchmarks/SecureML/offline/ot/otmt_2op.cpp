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

#define N 1024
//#define N 1280
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
SHOTIterated* ot;
SHOTIterated* ot2;

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


//SHOTExtension* ot;
//SHOTExtension* ot2;

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

int main(int argc, char** argv){
	
	srand (1);
	
	
	
	//setup connection
	int port, party;
	parse_party_and_port(argv, &party, &port);
	io = new NetIO(party==ALICE ? nullptr : "127.0.0.1", port);
	io->set_nodelay();

	ot = new SHOTIterated(io, party==ALICE);
	ot2 = new SHOTIterated(io, party!=ALICE);
	
	//ot = new SHOTExtension(io);
	//ot2 = new SHOTExtension(io);
	
	Mat a(N,D),b1(D,IT),c1(B,IT),b2(B,IT),c2(D,IT);
	
	vector<int> perm = random_perm();
	
	//for(int i=0;i<perm.size();i++)
	//	cout<<perm[i]<<",";
	//cout<<endl;
	
	
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
	Mat a_batch(B,D), a_batch_t(D,B), b1_IT(D,1),b2_IT(B,1);

	uint64 b_temp;
	
	uint64 r;
	
	Mat rsum(B,D);
	
	
	uint64 m0_temp,m1_temp;
	bool* select = new bool[(P+1)*B/4*D];
	block* m0 = new block[(P+1)*B/4*D], *m1 = new block[(P+1)*B/4*D], *mb = new block[(P+1)*B/4*D];
	
	int t2 = mytime.start();
	
	double ottime=0.0;
	

	
	
	for(int round=0;round<IT;round++){	
		
		int index = 0;
		//cout<<round<<endl;
		
		next_batch(a_batch, start_setup,perm,a);
		
		rsum.setZero();
		
		
		for(int j=0;j<D;j++){
			b_temp = b1(j,round);
			
			
			for(int k=0;k<P;k++){
				Mat m0_pack(P-k,1), m1_pack(P-k,1);
				m0_pack.setZero();
				m1_pack.setZero();
				
				
			
				for(int i=0;i<B/2-1;i++){

					r = (randomlong()<<k);
					rsum(i,j) += r;
					
					
					m0_temp = r;
					m1_temp = a_batch(i,j)*(uint64(1)<<k)+r;
					
					m0_temp>>=k;
					m1_temp>>=k;
					
					for(int t = 0;t<P-k;t++){
						m0_pack(t)+=(m0_temp&1);
						m0_pack(t)<<=1;
						m0_temp>>=1;
						m1_pack(t)+=(m1_temp&1);
						m1_pack(t)<<=1;
						m1_temp>>=1;
					
					}

				}
				
				r = (randomlong()<<k);
				rsum(B/2-1,j) += r;
				
				
				m0_temp = r;
				m1_temp = a_batch(B/2-1,j)*(uint64(1)<<k)+r;
				
				m0_temp>>=k;
				m1_temp>>=k;
				
				for(int t = 0;t<P-k;t++){
					m0_pack(t)+=(m0_temp&1);
					m0_temp>>=1;
					m1_pack(t)+=(m1_temp&1);
					m1_temp>>=1;
				
				}

				
				Mat m2_pack(P-k,1), m3_pack(P-k,1);
				m2_pack.setZero();
				m3_pack.setZero();
				
				for(int i=B/2;i<B-1;i++){

					r = (randomlong()<<k);
					rsum(i,j) += r;
					
					
					m0_temp = r;
					m1_temp = a_batch(i,j)*(uint64(1)<<k)+r;
					
					m0_temp>>=k;
					m1_temp>>=k;
					
					for(int t = 0;t<P-k;t++){
						m2_pack(t)+=(m0_temp&1);
						m2_pack(t)<<=1;
						m0_temp>>=1;
						m3_pack(t)+=(m1_temp&1);
						m3_pack(t)<<=1;
						m1_temp>>=1;
					
					}

				}
				
				
				
				r = (randomlong()<<k);
				rsum(B-1,j) += r;
				
				
				m0_temp = r;
				m1_temp = a_batch(B-1,j)*(uint64(1)<<k)+r;
				
				m0_temp>>=k;
				m1_temp>>=k;
				
				for(int t = 0;t<P-k;t++){
					m2_pack(t)+=(m0_temp&1);
					m0_temp>>=1;
					m3_pack(t)+=(m1_temp&1);
					m1_temp>>=1;
				
				}
				

				
				
				for(int t = 0; t<P-k;t++){
					m0[index+t] = makeBlock(m2_pack(t),m0_pack(t));
					m1[index+t] = makeBlock(m3_pack(t),m1_pack(t));
					select[index+t] = b_temp&1;
				}
				

				
				index += P-k;
				b_temp>>=1;
			}
			
			
			
			
		}
		
		block delta;
		PRG prg(fix_key);
		prg.random_block(&delta, 1);
		
		int t1 = mytime.start();
		
		if(party == ALICE){
		
			ot->send(m0,m1,(P+1)*B/4*D);
			ot2->recv(mb,select,(P+1)*B/4*D);
			
			//ot->send_cot(m0,delta,(P+1)*B/4*D);
			//ot2->recv_cot(mb,select,(P+1)*B/4*D);
			
		}
		else{
			
			ot->recv(mb,select,(P+1)*B/4*D);
			ot2->send(m0,m1,(P+1)*B/4*D);
			//ot2->recv_cot(mb,select,(P+1)*B/4*D);
			//ot->send_cot(m0,delta,(P+1)*B/4*D);
			
			
		}
		
		ottime+=mytime.end(t1);
		
		index = 0;
		Mat number_rec(B,1);
		number_rec.setZero();
		for(int j=0;j<D;j++){
			
			for(int k=0;k<P;k++){
				unsigned long mb0_temp,mb1_temp;
				Mat number(B,1);
				number.setZero();
				for(int t=P-k-1;t>0;t--){
					mb0_temp = ((unsigned long *)(&mb[index+t]))[0];
					mb1_temp = ((unsigned long *)(&mb[index+t]))[1];

					for(int i=B/2-1;i>=0;i--){
						number(i)+=(mb0_temp&1);
						number(i+B/2)+=(mb1_temp&1);
						number(i)<<=1;
						number(i+B/2)<<=1;
						mb0_temp>>=1;
						mb1_temp>>=1;
					}
				}
				mb0_temp = ((unsigned long *)(&mb[index]))[0];
				mb1_temp = ((unsigned long *)(&mb[index]))[1];

				for(int i=B/2-1;i>=0;i--){
					number(i)+=(mb0_temp&1);
					number(i+B/2)+=(mb1_temp&1);
					mb0_temp>>=1;
					mb1_temp>>=1;
				}
				
				for(int i=0;i<B;i++){
					number(i)<<=k;
					number_rec(i)+=number(i);
				}
				index+=P-k;
			}
			
		}
		
		c1.col(round) = number_rec-rsum.rowwise().sum()+a_batch*b1.col(round);
		start_setup+=B;
		
	}
	
	cout<<mytime.end(t2)<<"s!"<<endl;
	
	cout<<ottime<<"s      OT!"<<endl;
	
	return 0;
		
	Mat a_rec, b_rec, b2_rec, c_rec, c2_rec;
	Mat c1_IT(B,1),c2_IT(D,1);
	
	start_setup = 0;
	for(int round = 0;round<IT;round++){
		
		next_batch(a_batch,start_setup,perm,a);
		b1_IT = b1.col(round);
		b2_IT = b2.col(round);
		c1_IT = c1.col(round);
		c2_IT = c2.col(round);
		
		
		a_rec = reconstruct(a_batch,party);
		b_rec = reconstruct(b1_IT,party);
		//b2_rec = reconstruct(b2_IT,party);
		c_rec = reconstruct(c1_IT,party);
		//c2_rec = reconstruct(c2_IT,party);
		
		if(party==ALICE){
			//for(int i=0;i<a_rec.rows();i++){
			//	cout<<(a_rec.row(i)*b_rec).isApprox(c_rec.row(i))<<",";
			//}
			//cout<<endl;
			
			cout<<(a_rec*b_rec).isApprox(c_rec)<<"!!!   \n";
			
			//cout<<(a_rec.transpose()*b2_rec).isApprox(c2_rec)<<"!!!\n";
		}
		start_setup+=B;
	}
	
	
	delete ot;
	delete ot2;
	delete io;
	
	
}
