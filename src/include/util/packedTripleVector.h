#ifndef R_SA_LCP_PACKEDTRIPLEVECTOR_H
#define R_SA_LCP_PACKEDTRIPLEVECTOR_H
#include<sdsl/int_vector.hpp>

struct packedTripleVector {
    typedef uint64_t size_type;
    uint8_t a,b,c,width;
    sdsl::bit_vector data;

    packedTripleVector() = default;

    packedTripleVector(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t elements): a(a), b(b), c(c), width(a+b+c), data(width*elements,0) {
        if (a > 64 || b > 64 || c > 64) {
            std::cerr << "ELEMENTS MUST BE <= 64 BITS!" << std::endl;
            exit(1);
        }
    }

    template<unsigned el>
    inline std::pair<uint8_t, uint8_t> offW() const ;

    template<unsigned el>
    uint64_t get(const uint64_t ind) const {
        static_assert(el < 3, "The element accessed in a packed triple must be between 0 and 2 inclusive");
        uint8_t off, w;
        std::tie(off, w) = offW<el>();
        return data.get_int(ind*width + off, w);
    }

    template<unsigned el>
    void set(const uint64_t ind, const uint64_t val) {
        static_assert(el < 3, "The element accessed in a packed triple must be between 0 and 2 inclusive");
        uint8_t off, w;
        std::tie(off, w) = offW<el>();
        data.set_int(ind*width + off, val, w);
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(a, out, child, "a");
        bytes += sdsl::serialize(b, out, child, "b");
        bytes += sdsl::serialize(c, out, child, "c");
        bytes += sdsl::serialize(width, out, child, "width");
        bytes += sdsl::serialize(data, out, child, "data");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(a, in);
        sdsl::load(b, in);
        sdsl::load(c, in);
        sdsl::load(width, in);
        sdsl::load(data, in);
    }

    uint64_t size() const {
        return data.size()/width;
    }
};
            
template<>
inline std::pair<uint8_t, uint8_t> packedTripleVector::offW<0>()  const {
    return {0, a};
}

template<>
inline std::pair<uint8_t, uint8_t> packedTripleVector::offW<1>()  const {
    return {a, b};
}

template<>
inline std::pair<uint8_t, uint8_t> packedTripleVector::offW<2>()  const {
    return {a + b, c};
}

#endif //#ifndef R_SA_LCP_PACKEDTRIPLEVECTOR_H
