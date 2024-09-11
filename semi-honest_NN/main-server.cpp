#include "utils.hpp"
#include "NN_Server.hpp"

NetIO *io_server;
NetIO *io_client_alice;
NetIO *io_client_bob;
NetIO *io_client;

void setup(int party) {
    if (party == ALICE || party == BOB) {
        io_server = new NetIO(party==ALICE ? nullptr : "127.0.0.1", 1026);
        io_server->set_nodelay();
		setup_semi_honest(io_server, party);
		cerr << "semi honest set\n";
    }
	
    if (party == ALICE || party == EVE) {
        io_client_alice = new NetIO(party==EVE ? nullptr : "127.0.0.1", 1027);
        io_client_alice->set_nodelay();
    }
    if (party == BOB || party == EVE) {
        io_client_bob = new NetIO(party==EVE ? nullptr : "127.0.0.1", 1028);
        io_client_bob->set_nodelay();
    }
    
    if(party == ALICE) io_client = io_client_alice;
    if(party == BOB) io_client = io_client_bob;
    
}

void load_test_data(Mat& test_data, Mat& test_label) {
    ifstream infile2( "../data/mnist_test.csv" );
    int i=0;
    cout<<"load testing data.......\n";

    while(infile2) {
        string s;
        if (!getline(infile2,s)) break;
        istringstream ss(s);
        int temp;
        char c;
        
        //read label
        ss>>temp;
        ss>>c;
        test_label(i) = temp;//0~9
        
        //read data 
        for(int j=0; j<D; j++) {
            ss>>test_data(i,j);
            test_data(i, j) <<= L;
            ss>>c;
        }
        i++;
        if(i >= testN) break;
    }

    test_data.conservativeResize(i, D);
    test_label.conservativeResize(i, 1);

    infile2.close();
    return;
}



int main(int argc, char** argv) {
	srand((unsigned)time(NULL));
	int port, party;
    parse_party_and_port(argv, &party, &port);
    setup(party);
	
	Mat test_data(testN,D), test_label(testN,1);
	if(party == ALICE && EVALUATE) {
		load_test_data(test_data, test_label);
        cout << "test data loaded\n";
	}
	
	vector<int> perm(Ep*N,0);
	io_client->recv_data(&perm[0], sizeof(int) * perm.size());
	
	int hahaha1 = 666, hahaha2;
	if(party == ALICE) {
		io_server->send_data(&hahaha1, sizeof(int));
		io_server->flush();
		io_server->recv_data(&hahaha2, sizeof(int));
	}
	else if(party == BOB) {
		io_server->recv_data(&hahaha2, sizeof(int));
		io_server->send_data(&hahaha1, sizeof(int));
		io_server->flush();
	}
	
	NN_Server nn_server(party, io_server, io_client);
	nn_server.online(test_data, test_label, perm);
	
	return 0;
}
