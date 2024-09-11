#ifndef LAYER_HPP
#define LAYER_HPP

#include "utils.hpp"


static vec_t tenc[P];
	
class Layer {

public:

	Mat b_mask[L + 1];

	Mat W, W_grad; //(D_in, D_out)
	Mat bias, bias_grad;
	Mat z1, b1, tau1, select_comparison, z2, tau2;
	Mat in_x;
	int it;
	string id;
	long int lr;
	double COMM, COMM_client;
	
	NetIO * io_server;
	NetIO * io_client_alice;
	NetIO * io_client_bob;
	NetIO * io_client;
	int D_in, D_out;
	int party;
	
	
	PRF *prf_server, *prf_client, *prf_client_same, own_prf;
	
	double act_time, mul_time, off_time;
	Time_Calculator mytime;
	
	void set_io(NetIO *_io_server, NetIO *_io_client) {
		io_server = _io_server;
		io_client = _io_client;
	}
	
	
	void exchange_comparisons(Mat &temp_comparison) {
		int r = temp_comparison.rows(), c = temp_comparison.cols();
		int sz = (r * c + 63) >> 6;
		Mat temp(sz, 1);
		prf_server->cal(temp);
		for (int i = 0, cnt = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j, ++cnt) {
				int x = (cnt >> 6), y = (cnt & 63);
				temp_comparison(i, j) = (temp(x, 0) >> y) & 1;
			}
		}
	}

	Mat inner_product_input(Mat& z, Mat& x, Mat& tau) {
		Mat alpha = reconstruct(W - z, party, io_server);
		COMM += 2 * W.rows() * W.cols() * sizeof(unsigned long int);
		Mat result = (x * alpha) + tau;

		return result;
	}
	
	inline Mat inner_product(Mat x, Mat y, Mat e, Mat f, Mat& c) {
		
		Mat rs;
		if(party == ALICE) rs = x * f + e * (y - f) + c;
		else rs = x * f + e * y + c;
		return rs;
	}
	
	void update_W() {
		if(party == EVE) {
			cout << "clients update W???\n";
			assert(0);
		}
		//if(it % 30 == 0) lr = lr << 1;
		if(DEBUG_BACKWARD) {
			Mat rec_W_grad = reconstruct(W_grad, party, io_server);
			if(party == ALICE) debug_mat(rec_W_grad, id + "_W_grad.txt");
		}
		
		for (int i = 0; i < D_in; ++i) for (int j = 0; j < D_out; ++j) {
			W(i, j) = W(i, j) - ((long int)W_grad(i, j)) / (lr * B);
			if(it >= 50) {
				if(id == "L3") W(i, j) -= (long int)W(i, j) / 100;
				else W(i, j) -= (long int)W(i, j) / 400;
			}
		}
		for (int i = 0; i < D_out; ++i) {
			bias(0, i) = bias(0, i) - ((long int)bias_grad(0, i)) / (lr * B);
			if(it >= 50) {
				if(id == "L3") bias(0, i) -= (long int)bias(0, i) / 100;
				else bias(0, i) -= (long int)bias(0, i) / 400;
			}
		}
	}
	
/*------------------------------------------------------------------------*/
/*--------------------------activation function---------------------------*/
/*------------------------------------------------------------------------*/

	uint64 get_element(uint64 x, unsigned a, bool b, uint64 p, uint64 r) {
		uint64 modified;
		if ((unsigned) b == 0) {
		    if (a == 0) {
		        modified = (x + 1) << p;
		    } else {
		        modified = (1UL << (P - 1)) + r; // non-overlapping value
		    }
		} else {
		    if (a == 1) {
		        modified = x << p;
		    } else {
		        modified = (3UL << (P - 2)) + r; // represented as 10 in binary so non-overlapping value
		    }
		}
		return modified;
	}

	vec_t decompose_long(uint64 x, unsigned int length) {
		vec_t repr = new uint8_t[sizeof(uint64) * 2];
		for (size_t i = 0; i < sizeof(uint64) * 2; i++) {
		    repr[i] = (uint8_t) x;
		    x = x >> 8;
		}
		return repr;
	}

	void compute_encoding(bool deb, vec_t* vec, uint64 x, bool b) {
		unsigned a;
		uint64 modified;
		Mat r(P, 1);
		own_prf.cal(r);
		for (int i = 0; i < P; i++) {
		    a = x & 1;
		    modified = get_element(x, a, b, i, r(i, 0) >> 2);
		    if (deb)
		        cerr << i << ' ' << a << ' ' << b << ' ' << modified << ' ' << (r(i, 0) >> 2) << endl;
		    vec[i] = decompose_long(modified, P - i + 1);
		    x = x >> 1;
		}
	}

	void permute_encoding(vec_t* vec, PRP* prp) {
		for (int i = 0; i < P; i++) {
		    prp->permute_data(vec[i], sizeof(uint64) * 2);
		}
	}

	void cal_encodings(bool deb, uint64 share, bool comp, int party, PRP* prp, vec_t encoding[]) {
		
		if(DEBUG_ACT) {
			if(party == ALICE) {
				bool comp_;
				io_server->recv_data(&comp_, sizeof(bool));
				if(comp_ != comp) assert(0);
			}
			else {
				io_server->send_data(&comp, sizeof(bool));
				io_server->flush();
			}
		}
		
		if (party == ALICE) comp = !comp;

		compute_encoding(deb, &encoding[0], share, comp);
		permute_encoding(&encoding[0], prp);
	}

	static bool comp_(const int &x, const int &y) {
		for (int i = 0; i < 16; i++) {
		    if (tenc[x][i] < tenc[y][i]) return 1;
		    if (tenc[x][i] > tenc[y][i]) return 0;
		}
		return 0;
	}

	Mat cwise_product(Mat&x, Mat& z, Mat& y, int party) {
		
		Mat alpha = reconstruct(x - z, party, io_server);
		COMM += 2 * x.rows() * x.cols() * sizeof(unsigned long int);
		Mat result = alpha.cwiseProduct(y);

		return result;
	}
	
	Mat activation_function_server(Mat& d) {
		//cout << "in activation function ...\n";
		PRP* prp = new PRP();

		Mat d1(B, D_out);
		for (int i = 0; i < B; i++) for (int j = 0; j < D_out; ++j){
			// given d = d1 + d2 (shares, not the ones below) we want to find
			// d1 + d2 > 0 -> d1 > -d2
			if (party == ALICE) d1(i, j) = d(i, j);
			if (party == BOB) d1(i, j) = -d(i, j);
		}
		
		vec_t* encodings[B][D_out];
		for (int i = 0; i < B; i++) for (int j = 0; j < D_out; ++j) {
			encodings[i][j] = new vec_t[P];
			uint64 temp_x;
    		bool now, temp_comp;
    		
    		temp_x = d1(i, j);
			if(select_comparison(i, j) && party == BOB) temp_x = temp_x + 1;
			now = ((long int)d1(i, j) < 0);
			temp_comp = select_comparison(i, j) ^ now;
			if(now) temp_x = -temp_x;
    		cal_encodings(0, temp_x, temp_comp, party, prp, encodings[i][j]);
		}
		
		uint8_t* enc_t = new uint8_t[B * D_out * P * 16];
		int pid[P];
		int cnt = 0;
		for (int i = 0; i < B; i++) for (int j = 0; j < D_out; ++j) {
			for (int k = 0; k < P; k++) {
                pid[k] = k;
                tenc[k] = encodings[i][j][k];
            }
            sort(pid, pid + P, comp_);
			for (int k = 0; k < P; k++) for (int l = 0; l < 16; l++)
				enc_t[cnt++] = encodings[i][j][pid[k]][l];
		}
		
		io_client->send_data(&enc_t[0], sizeof(uint8_t) * B * D_out * P * 16);
		io_client->flush();
		
		COMM_client += sizeof(uint8_t) * B * D_out * P * 16;
		
		for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) for (int k = 0; k < P; k++)  {
			delete [] encodings[i][j][k];
		}
		for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) delete [] encodings[i][j];
		
		delete []enc_t;
		
		//cerr << "merging mats ...\n";
		COMM_client += merge_mat(b1, prf_client, io_client, party);
		COMM_client += merge_mat(tau1, prf_client, io_client, party);
		if(id != "L3") COMM_client += merge_mat(tau2, prf_client, io_client, party);
		
		Mat rs = cwise_product(d, z1, b1, party) + tau1;
		for (int i = 0; i < B; ++i) for (int j = 0; j < D_out; ++j) if(select_comparison(i, j)) 
			rs(i, j) = d(i, j) - rs(i, j);

		delete prp;
		return rs;
	}

	Mat sign_check(Mat& d, Mat& b1, Mat& s, Mat& b_mask) {
		//cout << "in activation function ...\n";
		PRP* prp = new PRP();

		int R = d.rows(), C = d.cols();

		Mat d1(R, C);
		for (int i = 0; i < R; i++)
			for (int j = 0; j < C; ++j){
			// given d = d1 + d2 (shares, not the ones below) we want to find
			// d1 + d2 >= 0 -> d1 + 1 > -d2
				if (party == ALICE) d1(i, j) = d(i, j) + 1;
				if (party == BOB) d1(i, j) = -d(i, j);
			}

		Mat select_comparison_(R, C);
		exchange_comparisons(select_comparison_);
		
		vec_t* encodings[R][C];
		for (int i = 0; i < R; i++)
			for (int j = 0; j < C; ++j) {
				encodings[i][j] = new vec_t[P];
				uint64 temp_x;
				bool now, temp_comp;
				
				temp_x = d1(i, j);
				if(select_comparison_(i, j) && party == BOB) temp_x = temp_x + 1;
				now = ((long int)d1(i, j) < 0);
				temp_comp = select_comparison_(i, j) ^ now;
				if(now) temp_x = -temp_x;
				cal_encodings(0, temp_x, temp_comp, party, prp, encodings[i][j]);
			}
		
		uint8_t* enc_t = new uint8_t[R * C * P * 16];
		int pid[P];
		int cnt = 0;
		for (int i = 0; i < R; i++) for (int j = 0; j < C; ++j) {
			for (int k = 0; k < P; k++) {
		    	pid[k] = k;
		        tenc[k] = encodings[i][j][k];
		    }
            sort(pid, pid + P, comp_);
			for (int k = 0; k < P; k++) for (int l = 0; l < 16; l++)
				enc_t[cnt++] = encodings[i][j][pid[k]][l];
		}
		//cerr << "sending encodings ...\n";
		io_client->send_data(&enc_t[0], sizeof(uint8_t) * R * C * P * 16);
		io_client->flush();
		
		COMM_client += sizeof(uint8_t) * R * C * P * 16;
				
		for (int i = 0; i < R; ++i) for (int j = 0; j < C; ++j) for (int k = 0; k < P; k++)  {
			delete [] encodings[i][j][k];
		}
		for (int i = 0; i < R; ++i) for (int j = 0; j < C; ++j) delete [] encodings[i][j];
		
		delete []enc_t;
		
		//cerr << "merging mats ...\n";
		Mat tau(R, C);
		COMM_client += merge_mat(b1, prf_client, io_client, party);
		COMM_client += merge_mat(tau, prf_client, io_client, party);
		for (int i = 0; i < R; ++i)
			for (int j = 0; j < C; ++j)
				if (select_comparison_(i, j)) {
					b1(i, j) = (party == ALICE) - b1(i, j);
					tau(i, j) = b_mask(i, j) - tau(i, j);
				}
		Mat rs = cwise_product(s, b_mask, b1, party) + tau;

		delete prp;
		return rs;
	}
	
	
	void print() {
		cout << "................... " << id << ":\n";
		cout << "offline: " << off_time << endl;
		cout << "multiplication: " << mul_time << endl;
		cout << "activation: " << act_time << endl;
		cout << "communication between servers:" << COMM / 1024 / 1024 << "MB\n";
		cout << "communication with clients:" << 2 * COMM_client / 1024 / 1024 << "MB\n";
	}
};

#endif
