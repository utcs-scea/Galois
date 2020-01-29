#ifndef _MODEL_H_
#define _MODEL_H_

#include <random>
#include "types.h"
#include "utils.h"
#include "lgraph.h"
#include "layers.h"
#include "math_functions.hpp"

#define NUM_LAYERS 2
std::string path = "/h2/xchen/datasets/Learning/"; // path to the input dataset

// N: number of vertices, D: feature vector dimentions, 
// E: number of distinct labels, i.e. number of vertex classes
// layer 1: features N x D, weights D x 16, out N x 16 (hidden1=16)
// layer 2: features N x 16, weights 16 x E, out N x E
class Model {
public:
	Model() {}

	// user-defined aggregate function
	virtual void aggregate(const VertexID src, const FV2D &features, FV &sum) {}
	
	// user-defined combine function
	virtual void combine(const FV2D ma, const FV2D mb, const FV &a, const FV &b, FV &out) {}
	
	void init() {
		read_graph(dataset, g); 
		n = g.size(); // N
		h_in.resize(n); // input embedding: N x D
		feature_dim = read_features(dataset, h_in);
		labels.resize(n, 0); // label for each vertex: N x 1
		num_classes = read_labels(dataset, labels);
		W.resize(NUM_LAYERS);
		Q.resize(NUM_LAYERS);
		init_matrix(feature_dim, hidden1, W[0]);
		init_matrix(hidden1, num_classes, W[1]);
		init_matrix(feature_dim, hidden1, Q[0]);
		init_matrix(hidden1, num_classes, Q[1]);

		train_mask.resize(n, 0);
		val_mask.resize(n, 0);
		//for (size_t i = 0; i < n; i ++) {
		//	train_mask[i] = 0;
		//	val_mask[i] = 0;
		//}
		set_masks(n, train_mask, val_mask);
		y_train.resize(n, 0);
		y_val.resize(n, 0);
		for (size_t i = 0; i < n; i ++) y_train[i] = (train_mask[i] == 1 ? labels[i] : -1);
		for (size_t i = 0; i < n; i ++) y_val[i] = (val_mask[i] == 1 ? labels[i] : -1);
		
		h_hidden1.resize(n); // hidden1 level embedding: N x 16
		h_out.resize(n); // output embedding: N x E
		h_softmax.resize(n); // normalized output embedding: N x E
		for (size_t i = 0; i < n; ++i) h_hidden1[i].resize(hidden1);
		for (size_t i = 0; i < n; ++i) h_out[i].resize(num_classes);
		for (size_t i = 0; i < n; ++i) h_softmax[i].resize(num_classes);
	}
	size_t get_nnodes() { return n; }
	size_t get_nedges() { return g.sizeEdges(); }
	size_t get_ft_dim() { return feature_dim; }
	size_t get_nclasses() { return num_classes; }
	size_t get_label(size_t i) { return labels[i]; }

	// forward pass
	void forward(LabelList labels, MaskList masks, AccT &loss, AccT &accuracy) {
		//size_t dim = h_in[0].size();
		//std::cout << "[debug] layer 0: input layer\n";
		galois::do_all(galois::iterate(g.begin(), g.end()), [&](const auto& src) {
			FV h_neighbors(feature_dim, 0); // used to gather neighbors' embeddings
			aggregate(src, h_in, h_neighbors);
			combine(W[0], Q[0], h_in[src], h_neighbors, h_hidden1[src]);
			relu(h_hidden1[src]);
		}, galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("Layer0"));

		//std::cout << "[debug] layer 1: hidden1 layer\n";
		galois::do_all(galois::iterate(g.begin(), g.end()), [&](const auto& src) {
			FV h_neighbors(hidden1, 0); // used to gather neighbors' embeddings
			aggregate(src, h_hidden1, h_neighbors);
			combine(W[1], Q[1], h_hidden1[src], h_neighbors, h_out[src]);
			relu(h_out[src]);
		}, galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("Layer1"));

		// TODO: need kernel fusion optimization
		// h_out is the output from output layer (num_examples x num_classes).
		// labels contains the ground truth labels for each example (num_examples x 1).
		// Note that labels is not one-hot encoded vector,
		// and it can be computed as y.argmax(axis=1) from one-hot encoded vector (y) of labels if required.
		std::vector<AccT> diffs(n, 0.0); // error for each vertex
		//std::vector<DataTy> y(num_classes); // ground truth
		//for (size_t src = 0; src < n; src++) {
		galois::do_all(galois::iterate(g.begin(), g.end()), [&](const auto& src) {
			if (labels[src] >= 0) { // masked
				softmax(h_out[src], h_softmax[src]); // normalize h_out using softmax
				// y  is a one hot encoded vector for the labels
				std::vector<AccT> y(num_classes, 0.0); // ground truth
				//for (size_t j = 0; j < num_classes; j ++) y[j] = 0.0; // ground truth
				y[labels[src]] = 1.0;
				diffs[src] = cross_entropy(y, h_softmax[src]);
			}
		}, galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("cross-entropy-back"));
		loss = masked_avg_loss(diffs, masks);
		//std::cout << "loss: " << loss << "\n";

		// comparing outputs (N x E) with the ground truth (labels)
		LabelList predictions(n);
		//auto num_classes = W[1][0].size();
		for (size_t i = 0; i < n; i ++) predictions[i] = argmax(num_classes, h_out[i]);
		accuracy = masked_accuracy(predictions, labels, masks);
	}

	// back propogation
	void backward(LabelList labels, MaskList masks) {
		FV2D in_diff(n);
		galois::do_all(galois::iterate(g.begin(), g.end()), [&](const auto& src) {
			in_diff[src].resize(num_classes);
			std::vector<AccT> y(num_classes, 0.0); // ground truth
			y[labels[src]] = 1.0;
			d_cross_entropy(y, h_softmax[src], in_diff[src]);
		}, galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("cross-entropy-back"));

		FV2D hidden1_diff(n);
		galois::do_all(galois::iterate(g.begin(), g.end()), [&](const auto& src) {
			FV temp_diff(num_classes);
			d_relu(in_diff[src], h_out[src], temp_diff);
			d_mvmul(temp_diff, hidden1_diff[src]);
		}, galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("layer1-back"));

		FV2D out_diff(n);
		galois::do_all(galois::iterate(g.begin(), g.end()), [&](const auto& src) {
			FV temp_diff(num_classes);
			d_relu(hidden1_diff[src], h_hidden1[src], temp_diff);
			d_mvmul(temp_diff, out_diff[src]);
		}, galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("layer0-back"));
	}

	// evaluate, i.e. inference or predict
	double evaluate(LabelList labels, MaskList masks, AccT &loss, AccT &acc) {
		Timer t_eval;
		t_eval.Start();
		auto num_classes = W[1][0].size();
		FV2D h_hidden1(n); // hidden1 level embedding
		FV2D h_out(n); // out level embedding
		for (size_t i = 0; i < n; ++i) h_hidden1[i].resize(hidden1);
		for (size_t i = 0; i < n; ++i) h_out[i].resize(num_classes);
		//std::cout << "loss: " << loss << "\n";
		forward(labels, masks, loss, acc);
		t_eval.Stop();
		return t_eval.Millisecs();
	}

	void train() {
		Timer t_epoch;
		// run epoches
		for (size_t i = 0; i < epochs; i++) {
			std::cout << "Epoch " << i << ": ";
			t_epoch.Start();
			// Construct feed dictionary

			// Training step
			AccT train_loss = 0.0, train_acc = 0.0;
			forward(y_train, train_mask, train_loss, train_acc);
			//backward(); // back propogation
			std::cout << " train_loss = " << train_loss << " train_acc = " << train_acc;

			// Validation
			AccT val_cost = 0.0, val_acc = 0.0;
			double eval_time = evaluate(y_val, val_mask, val_cost, val_acc);
			std::cout << " val_cost = " << val_cost << " val_acc = " << val_acc;

			t_epoch.Stop();
			std::cout << " time = " << t_epoch.Millisecs() << "ms. \n";
		}
	}

protected:
	size_t n; // N
	size_t feature_dim; // D
	size_t num_classes; // E
	Graph g; // the input graph
	FV2D h_in; // input_features: N x D
	FV2D h_hidden1; // hidden1 level embedding: N x 16
	FV2D h_out; // output embedding: N x E
	FV2D h_softmax; // output embedding: N x E
	std::vector<LabelT> labels; // labels for classification
	FV3D W; // parameters to learn, for vertex v, layer0: D x 16, layer1: 16 x E
	FV3D Q; // parameters to learn, for vertex u, i.e. v's neighbors, layer0: D x 16, layer1: 16 x E
	LabelList y_train, y_val; // labels for traning and validation
	MaskList train_mask, val_mask; // masks for traning and validation

	inline void init_matrix(size_t dim_x, size_t dim_y, FV2D &matrix) {
		// Glorot & Bengio (AISTATS 2010) init
		auto init_range = sqrt(6.0/(dim_x + dim_y));
		//std::cout << "Matrix init_range: (" << -init_range << ", " << init_range << ")\n";
		std::default_random_engine rng;
		std::uniform_real_distribution<FeatureT> dist(-init_range, init_range);
		matrix.resize(dim_x);
		for (size_t i = 0; i < dim_x; ++i) {
			matrix[i].resize(dim_y);
			for (size_t j = 0; j < dim_y; ++j)
				matrix[i][j] = dist(rng);
		}
		//for (size_t i = 0; i < 3; ++i)
		//	for (size_t j = 0; j < 3; ++j)
		//		std::cout << "matrix[" << i << "][" << j << "]: " << matrix[i][j] << std::endl;
	}

	inline void init_features(size_t dim, FV &x) {
		std::default_random_engine rng;
		std::uniform_real_distribution<FeatureT> dist(0, 0.1);
		for (size_t i = 0; i < dim; ++i)
			x[i] = dist(rng);
	}

	size_t read_labels(std::string dataset_str, LabelList &labels) {
		std::string filename = path + dataset_str + "-labels.txt";
		std::ifstream in;
		std::string line;
		in.open(filename, std::ios::in);
		size_t m, n;
		in >> m >> n >> std::ws;
		assert(m == labels.size()); // number of vertices
		std::cout << "label conuts: " << n << std::endl; // number of vertex classes
		IndexT v = 0;
		while (std::getline(in, line)) {
			std::istringstream label_stream(line);
			unsigned x;
			for (size_t idx = 0; idx < n; ++idx) {
				label_stream >> x;
				if (x != 0) {
					labels[v] = idx;
					break;
				}
			}
			v ++;
		}

		//for (size_t i = 0; i < 10; ++i)
		//	std::cout << "labels[" << i << "]: " << labels[i] << std::endl;
		return n;
	}

	size_t read_features(std::string dataset_str, FV2D &features) {
		std::string filename = path + dataset_str + ".ft";
		std::ifstream in;
		std::string line;
		in.open(filename, std::ios::in);
		size_t m, n;
		in >> m >> n >> std::ws;
		assert(m == features.size()); // m = number of vertices
		std::cout << "feature dimention: " << n << std::endl;
		for (size_t i = 0; i < m; ++i) {
			features[i].resize(n);
			for (size_t j = 0; j < n; ++j)
				features[i][j] = 0;
		}
		while (std::getline(in, line)) {
			std::istringstream edge_stream(line);
			IndexT u, v;
			FeatureT w;
			edge_stream >> u;
			edge_stream >> v;
			edge_stream >> w;
			features[u][v] = w;
		}
		/*
		for (size_t i = 0; i < 10; ++i) {
			for (size_t j = 0; j < n; ++j) {
				if (features[i][j] > 0)
					std::cout << "features[" << i << "][" << j << "]: " << features[i][j] << std::endl;
			}
		}
		//*/
		return n;
	}

	unsigned load_graph(Graph &graph, std::string filename, std::string filetype = "el") {
		LGraph lgraph;
		unsigned max_degree = 0;
		if (filetype == "el") {
			printf("Reading .el file: %s\n", filename.c_str());
			lgraph.read_edgelist(filename.c_str(), true); //symmetrize
			genGraph(lgraph, graph);
		} else if (filetype == "gr") {
			printf("Reading .gr file: %s\n", filename.c_str());
			galois::graphs::readGraph(graph, filename);
			galois::do_all(galois::iterate(graph.begin(), graph.end()), [&](const auto& vid) {
				graph.getData(vid) = 1;
				//for (auto e : graph.edges(n)) graph.getEdgeData(e) = 1;
			}, galois::chunk_size<256>(), galois::steal(), galois::loopname("assignVertexLabels"));
			std::vector<unsigned> degrees(graph.size());
			galois::do_all(galois::iterate(graph.begin(), graph.end()), [&](const auto& vid) {
				degrees[vid] = std::distance(graph.edge_begin(vid), graph.edge_end(vid));
			}, galois::loopname("computeMaxDegree"));
			max_degree = *(std::max_element(degrees.begin(), degrees.end()));
		} else { printf("Unkown file format\n"); exit(1); }
		if (filetype != "gr") {
			max_degree = lgraph.get_max_degree();
			lgraph.clean();
		}
		printf("max degree = %u\n", max_degree);
		return max_degree;
	}

	void genGraph(LGraph &lg, Graph &g) {
		g.allocateFrom(lg.num_vertices(), lg.num_edges());
		g.constructNodes();
		for (size_t i = 0; i < lg.num_vertices(); i++) {
			g.getData(i) = 1;
			auto row_begin = lg.get_offset(i);
			auto row_end = lg.get_offset(i+1);
			g.fixEndEdge(i, row_end);
			for (auto offset = row_begin; offset < row_end; offset ++)
				g.constructEdge(offset, lg.get_dest(offset), 0); // do not consider edge labels currently
		}
	}

	void read_graph(std::string dataset_str, Graph &g) {
		//printf("Start readGraph\n");
		galois::StatTimer Tread("GraphReadingTime");
		Tread.start();
		//std::string filename = dataset_str + ".gr";
		std::string filename = path + dataset_str + ".el";
		load_graph(g, filename);
		Tread.stop();
		//printf("Done readGraph\n");
		std::cout << "num_vertices " << g.size() << " num_edges " << g.sizeEdges() << "\n";
	}

	void set_masks(size_t n, MaskList &train_mask, MaskList &val_mask) {
		for (size_t i = 0; i < n; i++) {
			if (i < 120) train_mask[i] = 1; // [0, 120) train size = 120
			else if (i < 620) val_mask[i] = 1; // [120, 620) validation size = 500
			else ; // unlabeled vertices
		}
	}

};

#endif