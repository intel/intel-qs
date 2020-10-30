/**
 * @file Gates.hpp
 * @author Lee J. O'Riordan (lee.oriordan@ichec.ie)
 * @author Myles Doyle (myles.doyle@ichec.ie)
 * @brief Gates wrapper/storage class. Adapted from QNLP.
 * @version 0.2
 * @date 2020-06-01
 */

#ifndef GATECACHE
#define GATECACHE

#include <unordered_map>

#include <complex>
#include <vector>
#include <qureg.hpp>

#include <mat_ops.hpp>

#include "qureg.hpp"

/**
 * @brief Class to cache intermediate matrix values used within other parts of the computation. 
 * Heavily depended upon by NCU to store sqrt matrix values following Barenco et al. (1995) decomposition.
 * 
 * @tparam SimulatorType The simulator type with SimulatorGeneral as base class
 */
template <class Type>
class GateCache {
    private:
    //using GateType = TM2x2<Type>;
    std::size_t cache_depth;

    public:
    GateCache() : cache_depth(0) { };

    GateCache(std::size_t default_depth) : cache_depth(default_depth) {
        initCache(cache_depth);
    }

    ~GateCache(){ clearCache(); }

    //Take the 2x2 matrix type from the template SimulatorType
    using GateType = TM2x2<Type>;

    //Maintain a map for each gate label (X,Y,Z, etc.), and use vectors to store sqrt (indexed by 1/2^(i), and pairing matrix and adjoint)
    std::unordered_map<std::string, std::vector< std::pair<GateType, GateType> > > gateCacheMap;
    
    void clearCache(){
        gateCacheMap.clear();
        cache_depth = 0;
    }

    /**
     * @brief Initialise the gate cache with PauliX,Y,Z and H up to a given sqrt depth
     * 
     * @param sqrt_depth The depth to which calculate sqrt matrices and their respective adjoints
     */
    void initCache(const std::size_t sqrt_depth){
        // If we do not have a sufficient circuit depth, clear and rebuild up to given depth.

        if(cache_depth < sqrt_depth ){
            gateCacheMap.clear();
            cache_depth = 0;
        }

        if (gateCacheMap.empty()){
            gateCacheMap["X"] = std::vector< std::pair<GateType, GateType> > { std::make_pair( construct_pauli_x(), adjointMatrix( construct_pauli_x() ) ) };
            gateCacheMap["Y"] = std::vector< std::pair<GateType, GateType> > { std::make_pair( construct_pauli_y(), adjointMatrix( construct_pauli_y() ) ) };
            gateCacheMap["Z"] = std::vector< std::pair<GateType, GateType> > { std::make_pair( construct_pauli_z(), adjointMatrix( construct_pauli_z() ) ) };
            gateCacheMap["H"] = std::vector< std::pair<GateType, GateType> > { std::make_pair( construct_hadamard(), adjointMatrix( construct_hadamard() ) ) };

            for( std::size_t depth = 1; depth <= sqrt_depth; depth++ ){
                for( auto& kv : gateCacheMap ){
                    kv.second.reserve(sqrt_depth + 1);
                    auto m = matrixSqrt<GateType>(kv.second[depth-1].first);
                    kv.second.emplace(kv.second.begin() + depth, std::make_pair( m, adjointMatrix( m ) ) );
                }
            }
            cache_depth = sqrt_depth;
        }
    }

    /**
     * @brief Adds new gate to the cache up to a given sqrt depth
     * 
     * @param gateLabel Label of gate to index into map
     * @param gate Gate matrix
     * @param max_depth Depth of calculations for sqrt and associate adjoints
     */
    void addToCache(const std::string gateLabel, const GateType& gate, std::size_t max_depth){
        if(max_depth <= cache_depth && gateCacheMap.find(gateLabel) != gateCacheMap.end() ){
            return;
        }
        else if(max_depth > cache_depth){
            initCache(max_depth);
        }

        std::vector< std::pair<GateType, GateType> > v;

        v.reserve(max_depth + 1);
        v.push_back(std::make_pair( gate, adjointMatrix( gate ) ) );
        for( std::size_t depth = 1; depth <= max_depth; depth++ ){
            auto m = matrixSqrt<GateType>( v[depth-1].first );
            v.emplace(v.begin() + depth, std::make_pair( m, adjointMatrix( m ) ) );
        }
        gateCacheMap.emplace(std::make_pair(gateLabel, v) );
    }

    constexpr GateType construct_pauli_x(){
        GateType px;
        px(0, 0) = Type(0., 0.);
        px(0, 1) = Type(1., 0.);
        px(1, 0) = Type(1., 0.);
        px(1, 1) = Type(0., 0.);
        return px;
    }
    constexpr GateType construct_pauli_y(){
        GateType py;
        py(0, 0) = Type(0., 0.);
        py(0, 1) = Type(0., -1.);
        py(1, 0) = Type(0., 1.);
        py(1, 1) = Type(0., 0.);
        return py;
    }
    constexpr GateType construct_pauli_z(){
        GateType pz;
        pz(0, 0) = Type(1., 0.);
        pz(0, 1) = Type(0., 0.);
        pz(1, 0) = Type(0., 0.);
        pz(1, 1) = Type(-1., 0.);
        return pz;
    }
    constexpr GateType construct_hadamard(){
        GateType h;
        auto f = Type(1. / std::sqrt(2.), 0.0);
        h(0, 0) = h(0, 1) = h(1, 0) = f;
        h(1, 1) = -f;
        return h;
    }
};

template class GateCache<ComplexDP>;
template class GateCache<ComplexSP>;

#endif