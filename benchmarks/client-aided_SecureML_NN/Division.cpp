MT.cpp
	
	bool t;
    unsigned long rs, r1, r2;
    for (int j = 0; j < IT; j++) {
	    rs = rand();
	    r1 = r2 = 0;
	    for (int i = 63; i >= 0; --i) {
	        if ((rs >> i) & 1) {
	            t = rand() & 1;
	            if (t)
	                r2 += (1 << i);
	            else
	                r1 += (1 << i);
	        }
	    }  
	    F1 << rs << "," << r1 << "," << r2 << ",";
	    F2 << rs << "," << r1 << "," << r2 << ",";
	    F1 << endl;
	  	F2 << endl;
	}

main-server.cpp

	for (int i = 0; i < IT; i++) {
        string s;
        if (!getline(F,s)) {
            break;
        }
        istringstream ss(s);
        ss >> rs[i];
        ss >> temp;
        ss >> r1[i];
        ss >> temp;
        ss >> r2[i];
        ss >> temp;        
    }



	for (int i = 0; i < IT; ++i) {
        cout << "Iteration " << i << endl;

        int p = 1;
        Integer rs_ = Integer(64, rs[i], p);
        net->layers[LL]->set_r(rs_, party == ALICE ? r1[i] : r2[i]);


NeuralNetwork.hpp
	
	void computeGradient() {
	//	cout << "NN.computeGradient\n";

		int rows = B;
		int columns = LAST_LAYER_SIZE;
		int size = rows * columns;
		int index;

		Mat tmp_ = layers[LL]->y;
		Mat quotient(rows, columns);
		for (int i = 0; i < rows; ++i) {
			Batcher batcher1, batcher2;
			if (party == ALICE) {
				for (int j = 0; j < columns; ++j) {
					batcher1.add<Integer>(64, (unsigned long) tmp_(i, j));
	    			batcher2.add<Integer>(64, (unsigned long) 0);
				}
			}
			else {
				for (int j = 0; j < columns; ++j) {
					batcher1.add<Integer>(64, (unsigned long) 0);
	    			batcher2.add<Integer>(64, (unsigned long) tmp_(i, j));
				}
			}
			batcher1.make_semi_honest(ALICE);
	    	batcher2.make_semi_honest(BOB);
			Integer s = batcher1.next<Integer>();
			s = s + batcher2.next<Integer>();
			for (int j = 1; j < columns; ++j) {
				s = s + batcher1.next<Integer>();
				s = s + batcher2.next<Integer>();
			}

			for (int j = 0; j < columns; ++j) {
				Batcher batcher3, batcher4;
				if (party == ALICE) {
					batcher3.add<Integer>(64, (unsigned long) tmp_(i, j) * (1 << (L)));
	    			batcher4.add<Integer>(64, (unsigned long) 0);
				}
				else {
					batcher3.add<Integer>(64, (unsigned long) 0);
	    			batcher4.add<Integer>(64, (unsigned long) tmp_(i, j) * (1 << (L)));
				}
				batcher3.make_semi_honest(ALICE);
	    		batcher4.make_semi_honest(BOB);
	    		Integer temp1 = batcher3.next<Integer>();
				temp1 = temp1 + batcher4.next<Integer>();
				Integer ans = temp1 / s;
				int pp = 0;
				int p = 1;
				unsigned long x = ans.reveal<long long>(pp);
				Integer r = layers[LL]->rs;
				ans = ans ^ r;
				unsigned long x1 = ans.reveal<long long>(pp);
				unsigned long ra = layers[LL]->r;
				unsigned long ans_ = 0;
				for (int i = 63; i >= 0; --i) {
					unsigned long tmp0 = (party == ALICE ? (((unsigned long)x1 >> i) & 1) : 0) + ((ra >> i) & 1) - 2 * (((unsigned long)x1 >> i) & 1) * ((ra >> i) & 1);
					ans_ = (ans_ * 2) + tmp0;
				}
				quotient(i, j) = (unsigned long)ans_;
			}
		}

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < columns; ++j) {
				(*(layers[LL]->getGradient()))(i, j) = quotient(i, j) - outputData(i, j);
			}

		layers[LL]->fg = (*(layers[LL]->getGradient()));

		for (int i = LL; i > 0; --i) {
		//	cout << i << endl;
			layers[i]->computeGradient(*(layers[i-1]->getGradient()));
		}
	}
