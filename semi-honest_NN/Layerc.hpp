#ifndef LAYERC_HPP
#define LAYERC_HPP

#include "utils.hpp"


class Layerc {
public:
	string id;

	Mat b_mask[L + 1];

	Mat forward_z1, forward_z2;
	int D_in, D_out;
	int it;
	
	NetIO *io_alice, *io_bob;
	
	PRF *prf_alice, *prf_bob, *prf_alice_same, *prf_bob_same;
	
	void set_io(NetIO *_io_alice, NetIO *_io_bob) {
		io_alice = _io_alice;
		io_bob = _io_bob;
	}

	void set_id(string _id) {
		id = _id;
	}


/*------------------------------------------------------------------------*/
/*--------------------------activation function---------------------------*/
/*------------------------------------------------------------------------*/

	bool comp(const vec_t &repr1, const vec_t &repr2) {
		for (int i = 0; i < 16; i++) {
			if (repr1[i] < repr2[i]) return 1;
			if (repr1[i] > repr2[i]) return 0;
		}
		return 0;
	}

	bool check_reprs_equal(vec_t repr1, vec_t repr2, int length) {
		for (int i = 0; i < length; i++) {
			if (repr1[i] != repr2[i]) return false;
		}
		return true;
	}

	bool compute_vec_intersection(vec_t* vec1, vec_t* vec2) {
		unsigned int count = 0;
		int j = 0;
		for (int i = 0; i < P; i++){
			while (j < P && comp(vec2[j], vec1[i])) ++j;
			if (j >= P) break;
			if (check_reprs_equal(vec1[i], vec2[j], sizeof(unsigned long) * 2)) {
				count++;
				//cerr << i << ' ' << j << ' ' << vec1[i][0] << endl;
			}
		}
		assert(count == 0 || count == 1);
		return (bool)count;
	}
	
	void receive_encodings(Mat &counts) {
		//vector<uint8_t> enc0(B * D_out * P * 16, 0), enc1(B * D_out * P * 16, 0);
		uint8_t* enc0 = new uint8_t[B * D_out * P * 16];
		uint8_t* enc1 = new uint8_t[B * D_out * P * 16];
		//cout << "receiving encodings ...\n";
		io_alice->recv_data(&enc0[0], sizeof(uint8_t) * B * D_out * P * 16);
		io_bob->recv_data(&enc1[0], sizeof(uint8_t) * B * D_out * P * 16);
		//cout << "encodings received\n";
		
		vec_t alice_enc[P], bob_enc[P];
		for (int k = 0; k < P; ++k) {
			alice_enc[k] = new uint8_t[sizeof(unsigned long) * 2];
			bob_enc[k] = new uint8_t[sizeof(unsigned long) * 2];
		}
		int cnt = 0;
		for (int i = 0; i < B; i++) for (int j = 0; j < D_out; ++j) {
			for (int k = 0; k < P; k++) for (int l = 0; l < 16; l++) {
				alice_enc[k][l] = enc0[cnt];
				bob_enc[k][l] = enc1[cnt];
				++cnt;
			}	
			counts(i, j) = compute_vec_intersection(alice_enc, bob_enc);
		}
		delete[] enc0;
		delete[] enc1;
		for (int k = 0; k < P; k++) {
			delete[] alice_enc[k];
			delete[] bob_enc[k];
		}
	}

	void activation_function_client(Mat &z1, Mat &z2) {
		//cerr << "in activation function ...\n";
		Mat counts(B, D_out);
		Mat tau1(B, D_out), tau2(B, D_out);
		receive_encodings(counts);
		tau1 = counts.cwiseProduct(z1);
		if(id != "L3") tau2 = counts.cwiseProduct(z2);
		//cerr << "distributing mats ...\n";
		distribute_mat(counts, prf_alice, prf_bob, io_alice, io_bob);
		distribute_mat(tau1, prf_alice, prf_bob, io_alice, io_bob);
		if(id != "L3") distribute_mat(tau2, prf_alice, prf_bob, io_alice, io_bob);
	}

	void sign_check_client(Mat &z1) {
		Mat counts(B, D_out);
		Mat tau(B, D_out);
		receive_encodings(counts);
		tau = counts.cwiseProduct(z1);
		distribute_mat(counts, prf_alice, prf_bob, io_alice, io_bob);
		distribute_mat(tau, prf_alice, prf_bob, io_alice, io_bob);
	}

	
	void forward() {
		activation_function_client(forward_z1, forward_z2);
		++it;
	}

	void backward() {
		for (int i = L; i >= 0; i--)
			sign_check_client(b_mask[i]);
	}

};

#endif
