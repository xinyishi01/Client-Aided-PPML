#include "defines.hpp"

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

using namespace std;
using namespace Eigen;

class ReLU : public Layer {
public:
	Mat activation_function(Mat& inner, Mat& y, int party){
		Batcher batcher1, batcher2;
		
		if(party==ALICE){
			for(int i = 0; i < inner.rows(); ++i) {
				batcher1.add<Integer>(P, inner(i));
				batcher2.add<Integer>(P, (unsigned long)0);	
			}
		}
		else{
			for(int i = 0; i < inner.rows(); ++i) {
				batcher1.add<Integer>(P, (unsigned long)0);
				batcher2.add<Integer>(P, inner(i));	
			}
		}
		batcher1.make_semi_honest(ALICE);
		batcher2.make_semi_honest(BOB);
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s   1"<<endl;
		
		Bit* a = new Bit[inner.rows()];
		Bit* b = new Bit[inner.rows()];
		Bit* c = new Bit[inner.rows()];
		Bit* d = new Bit[inner.rows()];
		bool* lsb_c = new bool[inner.rows()];
		bool* lsb_d = new bool[inner.rows()];
		
		//Integer temp2(L+1,(int)pow(2,L),BOB);
		
		
		for(int i = 0; i < inner.rows(); ++i){
			//Integer temp = batcher1.next<Integer>() + batcher2.next<Integer>();
			Integer temp = batcher2.next<Integer>();
			Integer temp2 = batcher1.next<Integer>()+temp;
			b[i] = (temp2)[P-1];
			a[i] = (temp2 - Integer(P, 1ULL<<L, 0))[P-1] ;
			//a[i] = (batcher1.next<Integer>() + temp)[P-1];
			
			c[i] = a[i]&(!(b[i]));
			d[i] = !(a[i]);
			lsb_c[i] = getLSB(c[i]);
			lsb_d[i] = getLSB(d[i]);
			
		}
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s    2"<<endl;
		
		bool* select;
		select = new bool[inner.rows()*2];
		copy(lsb_c,lsb_c+inner.rows(),select);
		copy(lsb_d,lsb_d+inner.rows(),select+inner.rows());
		
		
		block* m0,*m1,*r;
		
		
		if(party==ALICE){
			
			m0 = new block[inner.rows()*2];
			m1 = new block[inner.rows()*2];
			r = new block[inner.rows()];
			
		}
		else{
			
			m0 = new block[inner.rows()];
			m1 = new block[inner.rows()];
			r = new block[inner.rows()*2];
		}
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 1"<<endl;
		
		uint64 temp_m0, temp_m1;
		Mat r1(inner.rows(), inner.cols()),r2,R(inner.rows(), inner.cols()),R2;
		
		for(int i=0;i<r1.rows();i++){
			r1(i) = randomlong();
		}
		
		if(party==ALICE){
			r2.resize(inner.rows(), inner.cols());
			
			for(int i=0;i<r2.rows();i++){
				r2(i) = randomlong();
			}
		}
		else{
			R2.resize(inner.rows(), inner.cols());
		}
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 2"<<endl;
		
		for(int i=0;i<inner.rows();++i){
			temp_m0 = (int)select[i]*inner(i)+r1(i);
			temp_m1 = (int)(!select[i])*inner(i)+r1(i);
			m0[i] = makeBlock((uint64)0,temp_m0);
			m1[i] = makeBlock((uint64)0,temp_m1);
		}
		
		
		if(party == ALICE){
			for(int i=inner.rows();i<inner.rows()*2;++i){
				temp_m0 = (int)select[i]*(uint64)(1<<L)+r2(i-inner.rows());
				temp_m1 = (int)(!select[i])*(uint64)(1<<L)+r2(i-inner.rows());
				m0[i] = makeBlock((uint64)(0),temp_m0);
				m1[i] = makeBlock((uint64)(0),temp_m1);
			}
		}
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 3"<<endl;
		
		Mat y_;
		
		if(party == ALICE){
		
			ot->send(m0,m1,inner.rows()*2);
			
			//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 4"<<endl;
			
			ot2->recv(r,select,inner.rows());
			
			for(int i=0;i<inner.rows();i++){	
				R(i) = ((unsigned long *)(&r[i]))[0];
			}
			
			y_ = R-r1-r2-y*(uint64)(1<<L);
			
		}
		else{
			ot->recv(r,select,inner.rows()*2);
			
			//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 4"<<endl;
			
			ot2->send(m0,m1,inner.rows());
			for(int i=0;i<inner.rows();i++){	
				R(i) = ((unsigned long *)(&r[i]))[0];
				R2(i) = ((unsigned long *)(&r[i+inner.rows()]))[0];
			}
			
			y_ = R+R2-r1-y*(uint64)(1<<L);
		}
		return y_;
	}

	void forward() {
		xa = reconstruct(x - a);
		U = MatMul(W, x, xa, b, c);
		y = activation_function(U);
	}
};