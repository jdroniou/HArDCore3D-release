#ifndef PARALLEL_FOR
#define PARALLEL_FOR

#include <thread>
#include <Eigen/Sparse>

namespace HArDCore3D
{
  
  /*!
   *	\addtogroup Common
   * @{
   */

  /// Struct to store the systems and vector
  template<typename MatrixType>
  struct SystemVectors
  {
    /// Constructor
    SystemVectors(std::vector<MatrixType> sys, std::vector<Eigen::VectorXd> vec):
      systems(sys),
      vectors(vec)
      {
      };
    
    std::vector<MatrixType> systems;
    std::vector<Eigen::VectorXd> vectors;
  };

  /// Function to distribute elements (considered as jobs) over threads. It returns a pair of vectors indicating the start and end element of each thread
  static std::pair<std::vector<int>, std::vector<int>> 
    distributeLoad(size_t nb_elements, unsigned nb_threads)
  { 
    // Vectors of start and end indices
    std::vector<int> start(nb_threads);
    std::vector<int> end(nb_threads);

    // Compute the batch size and the remainder
    unsigned batch_size = nb_elements / nb_threads;
    unsigned batch_remainder = nb_elements % nb_threads;

    // Distribute the remainder over the threads to get the start and end indices for each thread
    for (unsigned i = 0; i < nb_threads; ++i) {
      if (i < batch_remainder){
        start[i] = i * batch_size + i;
        end[i] = start[i] + batch_size + 1;
      }
      else{
        start[i] = i * batch_size + batch_remainder;
        end[i] = start[i] + batch_size;
      }
    }

    return std::make_pair(start, end);
  }

  /// Generic function to execute threaded processes
  static inline void parallel_for(unsigned nb_elements,
                                  std::function<void(size_t start, size_t end)> functor,
                                  bool use_threads = true)
  {
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    // Generate the start and end indices
    auto [start, end] = distributeLoad(nb_elements, nb_threads);

    std::vector<std::thread> my_threads(nb_threads);

    if (use_threads) {
      // Multithread execution
      for (unsigned i = 0; i < nb_threads; ++i) {
          my_threads[i] = std::thread(functor, start[i], end[i]);
      }
    } else {
      // Single thread execution (for easy debugging)
      for(unsigned i = 0; i < nb_threads; ++i) {
          functor(start[i], end[i]);
      }
    }

    // Wait for the other thread to finish their task
    if (use_threads) {
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
    }
  }

 /// Function to assemble global matrices from a procedure that compute local triplets
  static inline SystemVectors<Eigen::SparseMatrix<double>>
         parallel_assembly_system(
            size_t nb_elements,    //< nb of elements over which the threading will be done
            std::vector<std::pair<size_t, size_t>> size_systems,    //< sizes of each system to assemble
            std::vector<size_t> size_vectors,   //< sizes of each vector to assemble
            std::function<void(size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs)> batch_local_assembly, //< procedure to compute all the local contributions to the matrices (put in the triplets) and vectors between start and end points (e.g. elements indices) 
            bool use_threads = true   //< determine if threaded process is used or not
            )
  {
    // Matrices and vectors
    std::vector<Eigen::SparseMatrix<double>> systems;
    systems.reserve(size_systems.size());
    std::vector<Eigen::VectorXd> vectors;
    for (auto size : size_vectors){
      vectors.emplace_back(Eigen::VectorXd::Zero(size));
    }
    
    if (use_threads) {
      // Select the number of threads
      unsigned nb_threads_hint = std::thread::hardware_concurrency();
      unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

      // Generate the start and end indices
      auto [start, end] = distributeLoad(nb_elements, nb_threads);

      // Create vectors of triplets and vectors
      std::vector<std::vector<std::list<Eigen::Triplet<double> > > > triplets(nb_threads, std::vector<std::list<Eigen::Triplet<double> > >(size_systems.size()));
      std::vector<std::vector<Eigen::VectorXd>> vecs(nb_threads, vectors);

      // Assign a task to each thread
      std::vector<std::thread> my_threads(nb_threads);
      for (unsigned i = 0; i < nb_threads; ++i) {
          my_threads[i] = std::thread(batch_local_assembly, start[i], end[i], &triplets[i], &vecs[i]);
          // std::cout << "Thread " << i << ": " << start[i] << " " << end[i] << std::endl;
      }
    
      // Wait for the other threads to finish their task
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
    
      // Create systems from triplets
      for (size_t i = 0; i < size_systems.size(); i++){
        systems.emplace_back(Eigen::SparseMatrix<double>(size_systems[i].first, size_systems[i].second));
        size_t n_triplets = 0;
        for (auto triplets_thread : triplets) {
          n_triplets += triplets_thread[i].size();
        }
        std::vector<Eigen::Triplet<double>> all_triplets(n_triplets);
        auto triplet_index = all_triplets.begin();
        for (auto triplets_thread : triplets) {
          triplet_index = std::copy(triplets_thread[i].begin(), triplets_thread[i].end(), triplet_index);
        }
        systems[i].setFromTriplets(all_triplets.begin(), all_triplets.end());
      }
      
      for (size_t i = 0; i < size_vectors.size(); i++){
        for (auto vec_thread : vecs){
          vectors[i] += vec_thread[i];
        }
      }

    } else {
      std::vector<std::list<Eigen::Triplet<double> > > triplets(size_systems.size());
      batch_local_assembly(0, nb_elements, &triplets, &vectors);
      for (size_t i = 0; i < size_systems.size(); i++){
        systems.emplace_back(Eigen::SparseMatrix<double>(size_systems[i].first, size_systems[i].second));
        systems[i].setFromTriplets(triplets[i].begin(), triplets[i].end());
      }
    }
    return SystemVectors(systems, vectors);  
  }
  
  /// Function to assemble a global matrix and right-hand side from a procedure that compute local triplets and rhs contributions (a wrapper for the more general function that can assemble several matrices and vectors)
  static inline std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
         parallel_assembly_system(
            size_t nb_elements, //< nb of elements over which the threading will be done
            size_t size_system, //< size of the system to assemble
            std::function<void(size_t start, size_t end, std::list<Eigen::Triplet<double>> * triplets, Eigen::VectorXd * rhs)> batch_local_assembly, //< procedure to compute all the local contributions to the matrix (put in the triplets) and rhs (put in the vector) between start and end points (e.g. elements indices) 
            bool use_threads = true //< determine if threaded process is used or not
        )
  {
    std::function<void(size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs)> 
    new_batch_local_assembly = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs)-> void
    {batch_local_assembly(start, end, &(*triplets)[0], &(*vecs)[0]);};
    auto [systems, vectors] = parallel_assembly_system(nb_elements, {std::make_pair(size_system, size_system)}, {size_system}, new_batch_local_assembly, use_threads);
    return std::make_pair(systems[0], vectors[0]); 
  }
  
    
  /// Function to assemble two global matrices and vectors (such as: system and static condensation operator, or system and matrix for BC) from a procedure that compute local triplets and rhs contributions  (a wrapper for the more general function that can assemble several matrices and vectors)
  static inline std::tuple<Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::SparseMatrix<double>, Eigen::VectorXd>
         parallel_assembly_system(
            size_t nb_elements, //< nb of elements over which the threading will be done
            size_t size_system1, //< size of the system 1 to assemble (must be square, corresponds to first matrix and vector=rhs)
            std::pair<size_t, size_t> size_Mat2, //< size of the second matrix to assemble (can be rectangular)
            size_t size_b2, //< size of second vector
            std::function<void(size_t start, size_t end, std::list<Eigen::Triplet<double>> * triplets1, Eigen::VectorXd * vec1, std::list<Eigen::Triplet<double>> * triplets2, Eigen::VectorXd * vec2)> batch_local_assembly, //< procedure to compute all the local contributions to the matrices (put in the triplets) and vectors between start and end points (e.g. elements indices) 
            bool use_threads = true //< determine if threaded process is used or not
        )
  {
    std::function<void(size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs)> 
    new_batch_local_assembly = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs)-> void
    {batch_local_assembly(start, end, &(*triplets)[0], &(*vecs)[0], &(*triplets)[1], &(*vecs)[1]);};
    auto [systems, vectors] = parallel_assembly_system(nb_elements, {std::make_pair(size_system1, size_system1), size_Mat2}, {size_system1, size_b2}, new_batch_local_assembly, use_threads);
    return std::make_tuple(systems[0], vectors[0], systems[1], vectors[1]);
  }
  //@}
  
} // end of namespace HArDCore3D
#endif
