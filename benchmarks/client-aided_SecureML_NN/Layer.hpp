#include "defines.hpp"

using namespace std;
using namespace Eigen;

class Layer {
public:

	int t_mat1, t1;

	Time_Calculater mytime;

	NetIO *io;
	SHOTIterated *ot;
	SHOTIterated *ot2;

	Mat a1, b1[B], c1; // a1:B*d(i-1) b1:d(i-1)*d(i) c1:B*d(i)

	Mat a2, b2, c2; // a2,b2,c2:B*d(i)

	Mat a3, b3[B], c3; // a3:B*d(i) b3:d(i)*d(i-1) c3:B*d(i-1)

	Mat a4, b4, c4; // a4:d(i-1)*B b4:B*d(i) c4:d(i-1)*d(i)

	Mat x, xa, U, y; // x:B*d(i-1) U = x * W + bias(B*d(i)) y:B*d(i)

	Mat G, derivative_activation_function; //G,...:B*d(i)

	unsigned long r;
	Integer rs;

	int inputDim, outputDim;
	int party;
	int is_first;

	long int lr;

	Mat W, bias; // W:d(i-1)*d(i) bias:1*d(i)


	Layer(int _party) {
	//	initialize(_party);
		party = _party;
		is_first = 0;
		t_mat1 = 0;
	}

	myType myrand(int x) { // -len/10000 ~ len/10000
	//	return (5 * (x % 21)) * (1 << L) / 10000 - 50 * (1 << L) / 10000;
		const unsigned long O = 50, H = 50;
		unsigned long int rd = rand() % (O + H);
		rd = rd * (1 << L) / 10000;
		return rd - H * (1 << L) / 10000;
	}

	void initialize(myType _lr) {
		lr = _lr;
		W.resize(inputDim, outputDim);
		for (int i = 0; i < inputDim; ++i)
			for (int j = 0; j < outputDim; ++j)
				W(i,j) = myrand(i+j);
		bias.resize(1, outputDim);
		for (int i = 0; i < outputDim; ++i)
			bias(0,i) = 0;
		G.resize(B, outputDim);
		y.resize(B, outputDim);
		derivative_activation_function.resize(B, outputDim);	
	}
/*
	void set_MT(Mat _a1, Mat _b1, Mat _c1, Mat _a2, Mat _b2, Mat _c2, Mat _a3, Mat _b3, Mat _c3, Mat _a4, Mat _b4, Mat _c4) {
		a1 = _a1;
		b1 = _b1;
		c1 = _c1;
		a2 = _a2;
		b2 = _b2;
		c2 = _c2;
		a3 = _a3;
		b3 = _b3;
		c3 = _c3;
		a4 = _a4;
		b4 = _b4;
		c4 = _c4;
	}*/

	void set_r(Integer _rs, unsigned long _r) {
		rs = _rs;
		r = _r;
	}

	void set_io(NetIO *_io, SHOTIterated *_ot, SHOTIterated *_ot2) {
		io = _io;
		ot = _ot;
		ot2 = _ot2;
	}

	void set_xa(Mat _xa) {
		xa = _xa;
	}

	myType randomlong(){
		myType rand1 = abs(rand());
	    myType rand2 = abs(rand());
	    rand1 = rand1 << (sizeof(int)*8);   
	    myType randULL = (rand1 | rand2);   
	    return randULL;
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

	Mat reconstruct(Mat A) {
		vector<unsigned long int> A_temp(A.cols()*A.rows());
		
		if(party==ALICE){
			io->recv_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
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
			io->send_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
			io->flush();
			
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
			io->send_data(&A_temp[0],sizeof(unsigned long int)*A_temp.size());
			io->flush();
			
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
	}

	Mat spe_product(Mat a, Mat *b) {
		Mat c(B, b[0].cols());
		for (int i = 0; i < B; i++) {
			c.row(i) = a.row(i) * b[i];
		}
		return c;
	}

	Mat compute_inner0(Mat& W0,Mat& x0, Mat& e, Mat& b0,Mat& c0, Mat& Wb1) {
	    Mat f,d;
	    f = W0-b0+Wb1;
	    d = x0*f+e*W0+c0;
	    return d;
	}

	Mat spe_compute_inner0(Mat& W0,Mat& x0, Mat& e, Mat *b0,Mat& c0, Mat *Wb1) {
		Mat f[B],d;
		for (int i = 0; i < B; i++)
			f[i] = Wb1[i] - b0[i] + W0;
		d = spe_product(x0, f)+e*W0+c0;
		return d;
	}

	Mat compute_inner1(Mat& W1,Mat& x1, Mat& e, Mat& b1,Mat& c1, Mat& Wb0) {
	    Mat f,d;
	    f = W1-b1+Wb0;
	    d = x1*f+e*(W1-f)+c1;
	    return d;
	}

	Mat spe_compute_inner1(Mat& W1,Mat& x1, Mat& e, Mat *b1,Mat& c1, Mat *Wb0) {
	    Mat f[B], f_[B], d;
	    for (int i = 0; i < B; i++)
			f[i] = f_[i] = Wb0[i] - b1[i];
		for (int i = 0; i < B; i++) {
			f[i] = f[i] + W1;
			f_[i] = W1 - f[i];
		}
		d = spe_product(x1, f) + spe_product(e, f_) + c1;
		return d;
	}

	Mat MatMul(Mat W, Mat x, Mat a_IT1, Mat b_IT1, Mat c_IT1) {
    	Mat xa(x.rows(), x.cols()), W0b(W.rows(), W.cols()), W1b(W.rows(), W.cols());

    	xa = reconstruct(x - a_IT1);

    	W0b = W1b = W - b_IT1;
    	if (party == ALICE) {
    		send_mat(W0b, io);
    		receive_mat(W1b, io);
    	}
    	else {
    		receive_mat(W1b, io);
    		send_mat(W0b, io);
    	}

		Mat inner;

	//	cerr << xa(0, 0) << ' ' << W1b(0, 0) << endl;

		if(party == ALICE) {
			inner = compute_inner0(W, x, xa, b_IT1, c_IT1, W1b);
		}
		else {
			inner = compute_inner1(W, x, xa, b_IT1, c_IT1, W1b);
		}
	//	cout << "after com_inner\n";
		return inner;
	}

	Mat spe_MatMul(Mat W, Mat x, Mat a_IT1, Mat *b_IT1, Mat c_IT1) {
    	Mat xa(x.rows(), x.cols()), W0b[B], W1b[B];

    	xa = reconstruct(x - a_IT1);

    //	cerr << "aft xa\n";

    	for (int i = 0; i < B; i++)
    		W0b[i] = W1b[i] = W - b_IT1[i];
    	if (party == ALICE) {
    		for (int i = 0; i < B; i++) {
    		//	cerr << i << endl;
    			send_mat(W0b[i], io);
    		}
    		for (int i = 0; i < B; i++) {
    		//	cerr << i << endl;
    			receive_mat(W1b[i], io);
    		}
    	}
    	else {
    		for (int i = 0; i < B; i++) {
    		//	cerr << i << endl;
    			receive_mat(W1b[i], io);
    		}
    		for (int i = 0; i < B; i++) {
    		//	cerr << i << endl;
    			send_mat(W0b[i], io);
    		}
    	}

    //	cerr << "aft w-b\n";

		Mat inner;

	//	cerr << xa(0, 0) << ' ' << W1b(0, 0) << endl;

		if(party == ALICE) {
			inner = spe_compute_inner0(W, x, xa, b_IT1, c_IT1, W1b);
		}
		else {
			inner = spe_compute_inner1(W, x, xa, b_IT1, c_IT1, W1b);
		}
	//	cout << "after com_inner\n";
		return inner;
	}

	Mat MatMul_(Mat W, Mat x, Mat xa, Mat b_IT1, Mat c_IT1) {
    	Mat W0b(W.rows(), W.cols()), W1b(W.rows(), W.cols());

		if(party == ALICE){
			W0b = W-b_IT1;
		}
		else{
			W1b = W-b_IT1;
		}
		
		vector<unsigned long int> temp0(W0b.cols()*W0b.rows()),temp1(W1b.cols()*W1b.rows());
	
		if(party == ALICE) {
			for(int j=0;j<W0b.rows();j++) {
				for(int k=0;k<W0b.cols();k++)
					temp0[j*W0b.cols()+k] = W0b(j,k);
			}
			
			io->send_data(&temp0[0],sizeof(unsigned long int)*temp0.size());
			io->flush();
			io->recv_data(&temp1[0],sizeof(unsigned long int)*temp1.size());
			
			for(int j=0;j<W1b.rows();j++) {
				for(int k=0;k<W1b.cols();k++)
					W1b(j,k) = temp1[j*W1b.cols()+k];
			}
		}
		else {
			for(int j=0;j<W1b.rows();j++) {
				for(int k=0;k<W1b.cols();k++)
					temp1[j*W1b.cols()+k] = W1b(j,k);
			}
			
			io->recv_data(&temp0[0],sizeof(unsigned long int)*temp0.size());
			io->send_data(&temp1[0],sizeof(unsigned long int)*temp1.size());
			io->flush();
		
			for(int j=0;j<W0b.rows();j++){
				for(int k=0;k<W0b.cols();k++)
					W0b(j,k) = temp0[j*W0b.cols()+k];
			}
		}

		Mat inner;

		if(party == ALICE) {
			inner = compute_inner0(W, x, xa, b_IT1, c_IT1, W1b);
		}
		else {
			inner = compute_inner1(W, x, xa, b_IT1, c_IT1, W0b);
		}
		return inner;
	}

	Mat activation_function(Mat& inner, int p){
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
		
	//	Bit* a = new Bit[inner.rows()];
		Bit* b = new Bit[inner.rows()];
		Bit* c = new Bit[inner.rows()];
	//	Bit* d = new Bit[inner.rows()];
		bool* lsb_c = new bool[inner.rows()];
	//	bool* lsb_d = new bool[inner.rows()];
		
		//Integer temp2(L+1,(int)pow(2,L),BOB);
		
		for(int i = 0; i < inner.rows(); ++i){
			//Integer temp = batcher1.next<Integer>() + batcher2.next<Integer>();
			Integer temp = batcher2.next<Integer>();
			Integer temp2 = batcher1.next<Integer>()+temp;
			b[i] = (temp2)[P-1];
		//	a[i] = (temp2 - Integer(P, 1ULL<< (L), 0))[P-1];

			c[i] = !(b[i]);
		//	d[i] = 0;
			lsb_c[i] = getLSB(c[i]);
		//	lsb_d[i] = getLSB(d[i]);	
		}

		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s    2"<<endl;
		
		bool* select;
		select = new bool[inner.rows() * 2];
		copy(lsb_c,lsb_c+inner.rows(),select);
		copy(lsb_c,lsb_c+inner.rows(),select+inner.rows());
		
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
		
		myType temp_m0, temp_m1;
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
			m0[i] = makeBlock((myType)0,temp_m0);
			m1[i] = makeBlock((myType)0,temp_m1);
		}
		
		if(party == ALICE){
			for(int i=inner.rows();i<inner.rows()*2;++i){
				temp_m0 = (int)select[i - inner.rows()]*(myType)(1<< (L))+r2(i-inner.rows());
				temp_m1 = (int)(!select[i - inner.rows()])*(myType)(1<< (L))+r2(i-inner.rows());
				m0[i] = makeBlock((myType)(0),temp_m0);
				m1[i] = makeBlock((myType)(0),temp_m1);
			}
		}
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 3"<<endl;
		
		Mat y_;
		Mat y__;
		
		if(party == ALICE){
		
			ot->send(m0,m1,inner.rows()*2);
			
			//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 4"<<endl;
			
			ot2->recv(r,select,inner.rows());
			
			for(int i=0;i<inner.rows();i++){	
				R(i) = ((unsigned long *)(&r[i]))[0];
			}
			
			y_ = R - r1;

			y__ = -r2;
			
		}
		else{
			ot->recv(r,select,inner.rows()*2);
			
			//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s 4"<<endl;
			
			ot2->send(m0,m1,inner.rows());

			for(int i=0;i<inner.rows();i++){	
				R(i) = ((unsigned long *)(&r[i]))[0];
				R2(i) = ((unsigned long *)(&r[i+inner.rows()]))[0];
			}
			
			y_ = R - r1;

			y__ = R2;
		}

	//	cerr << y__.rows() << ' ' << y__.cols() << endl;
	//	cerr << derivative_activation_function.rows() << ' ' << 1 << endl;

		derivative_activation_function.col(p) = y__;

		return y_;
	}

	Mat get_derivative(Mat& inner){
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
		
		for(int i = 0; i < inner.rows(); ++i){
			//Integer temp = batcher1.next<Integer>() + batcher2.next<Integer>();
			Integer temp = batcher2.next<Integer>();
			Integer temp2 = batcher1.next<Integer>()+temp;
			b[i] = (temp2)[P-1];
			a[i] = (temp2 - Integer(P, 1ULL<< (L), 0))[P-1];
			
			c[i] = 0;
			d[i] = !(b[i]);
			lsb_c[i] = getLSB(c[i]);
			lsb_d[i] = getLSB(d[i]);
			
		}
		
		//cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s2"<<endl;
		
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
		
		myType temp_m0, temp_m1;
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
		
		for(int i=0;i<inner.rows();++i){
			temp_m0 = (int)select[i]+r1(i);
			temp_m1 = (int)(!select[i])+r1(i);
			m0[i] = makeBlock((myType)0,temp_m0);
			m1[i] = makeBlock((myType)0,temp_m1);
		}
		
		
		if(party == ALICE){
			for(int i=inner.rows();i<inner.rows()*2;++i){
				temp_m0 = (int)select[i]*(myType)(1<< (L))+r2(i-inner.rows());
				temp_m1 = (int)(!select[i])*(myType)(1<< (L))+r2(i-inner.rows());
				m0[i] = makeBlock((myType)(0),temp_m0);
				m1[i] = makeBlock((myType)(0),temp_m1);
			}
		}
		
		Mat y_;
		
		if(party == ALICE){
		
			ot->send(m0,m1,inner.rows()*2);
			
			ot2->recv(r,select,inner.rows());
			
			for(int i=0;i<inner.rows();i++){	
				R(i) = ((unsigned long *)(&r[i]))[0];
			}
			
			y_ = R-r1-r2;
			
		}
		else{
			ot->recv(r,select,inner.rows()*2);
			
			ot2->send(m0,m1,inner.rows());
			for(int i=0;i<inner.rows();i++){	
				R(i) = ((unsigned long *)(&r[i]))[0];
				R2(i) = ((unsigned long *)(&r[i+inner.rows()]))[0];
			}
			
			y_ = R+R2-r1;
		}
		return y_;
	}

	void forward(const Mat& x_input) {
		x = x_input;
	//	cerr << "bef Matmul\n";

		if (EVAT)
			t1 = mytime.start();
/*
		if (is_first)
			U = MatMul_(W, x, xa, b1, c1);
		else	*/
		U = spe_MatMul(W, x, a1, b1, c1);

		if (EVAT)
			cerr << "MatMul1: " << mytime.end(t1) << endl;

	//	cout << "after 1 spe_MatMul\n";
		trunc(U);

		for (int i = 0; i < B; ++i)
			U.row(i) = U.row(i) + bias;

		if (EVAT)
			t1 = mytime.start();

		for (int i = 0; i < outputDim; ++i) {
			Mat tmp = U.col(i);
		//	cout << tmp.rows() << ' ' << tmp.cols() << endl;
			if (ACT)
				y.col(i) = activation_function(tmp, i);
		}

		if (EVAT)
			cerr << "Act1: " << mytime.end(t1) << endl;
	}

	Mat MatMul_cwise(Mat Wb, Mat xa, Mat W, Mat x, Mat a_IT1, Mat b_IT1, Mat c_IT1) {
		Mat inner;

	//	cerr << xa(0, 0) << ' ' << Wb(0, 0) << endl;

		if(party == ALICE) {
			inner = compute_inner0(W, x, xa, b_IT1, c_IT1, Wb);
		}
		else {
			inner = compute_inner1(W, x, xa, b_IT1, c_IT1, Wb);
		}
	//	cout << "after com_inner\n";
		return inner;
	}

	void computeGradient(Mat& prevGradient) {

		if (EVAT)
			t1 = mytime.start();

		int rows = B;
		int columns = inputDim;
		int common_dim = outputDim;
		int tempSize = rows*common_dim;

		Mat tmp(rows, common_dim);
		Mat tG(1, 1), td(1, 1), ta2(1, 1), tb2(1, 1), tc2(1, 1);
/*
		Mat derivative_activation_function_;
		derivative_activation_function_.resize(B, outputDim);
		for (int i = 0; i < outputDim; ++i) {
			Mat tmp_ = U.col(i);
		//	cout << tmp.rows() << ' ' << tmp.cols() << endl;
			if (ACT)
				derivative_activation_function_.col(i) = get_derivative(tmp_);
			
			if (reconstruct(derivative_activation_function_.col(i)) != reconstruct(derivative_activation_function.col(i)))
				cerr << i << ' ' << "!!!!!!!!!!!!\n";
		}*/

		if (outputDim != 10) {
			tmp = G;
			Mat Wb = G - b2, Wb_, xa = derivative_activation_function - a2;
			xa = reconstruct(xa);
			Wb_ = Wb;
			if (party == ALICE) {
				receive_mat(Wb_, io);
				send_mat(Wb, io);
			}
			else {
				send_mat(Wb, io);
				receive_mat(Wb_, io);
			}
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < common_dim; ++j) {
				//	cout << i << ' ' << j << endl;
					Matrix<myType, 1, 1> tG{tmp(i, j)};
					Matrix<myType, 1, 1> td{derivative_activation_function(i, j)};
					Matrix<myType, 1, 1> ta2{a2(i, j)};
					Matrix<myType, 1, 1> tb2{b2(i, j)};
					Matrix<myType, 1, 1> tc2{c2(i, j)};
					Matrix<myType, 1, 1> tWb_{Wb_(i, j)};
					Matrix<myType, 1, 1> txa{xa(i, j)};
					G(i, j) = MatMul_cwise(tWb_, txa, tG, td, ta2, tb2, tc2)(0, 0);
				//	cerr << i << ' ' << j << ' ' << G(i, j) << ' ';
				//	G(i, j) = MatMul(tG, td, ta2, tb2, tc2)(0, 0);
				//	cerr << "!!!!!!!!!!\n";
				//	return;
				//	cerr << G(i, j) << endl;
				}
			trunc(G);
		}

		if (EVAT)
			cerr << "cwiseMat: " << mytime.end(t1) << endl;

	//	cout << "-------before MatMul\n";

		if (EVAT)
			t1 = mytime.start();

		if (inputDim != LAYER0)
			prevGradient = spe_MatMul(W.transpose(), G, a3, b3, c3);

		if (EVAT)
			cerr << "matmul3: " << mytime.end(t1) << endl;

		trunc(prevGradient);

	//	cout << "after MatMul\n";

	}

	void updateEquations(const Mat& prevActivations) {
		//Update Bias
		myType sum;

		for (int i = 0; i < outputDim; ++i){
			sum = 0;
			for (int j = 0; j < B; ++j)	
				sum += G(j, i);
		//	cout << i << endl;
			bias(0, i) -= (long int)sum / ((long int)lr*B) + (long int)bias(0, i) / K;
		}

	//	cout << "after bias\n";

		//Update Weights
		int rows = inputDim;
		int columns = outputDim;
		int common_dim = B;

		Mat delta_W(rows, columns);
		
		delta_W = MatMul(G, prevActivations.transpose(), a4, b4, c4);

		trunc(delta_W);

	//	cout << "after MatMul\n";

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < columns; ++j)
				delta_W(i, j) = (long int)delta_W(i, j) / ((long int)lr*B);

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < columns; ++j)
				W(i, j) -= delta_W(i, j) + (long int)W(i, j) / K;
	}

	//getters

	Mat* getActivation() {
		return &y;
	}

	Mat* getGradient() {
		return &G;
	}
};
