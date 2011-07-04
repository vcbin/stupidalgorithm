#ifndef STUPIDALGO_ALL_ALG_TYPEDEF
#define STUPIDALGO_ALL_ALG_TYPEDEF

enum alg_type
{
	alg_pso=1,alg_mpso,alg_arpso,alg_dpso,alg_dpso_m,alg_psobc,
	alg_de,alg_bbde,alg_sde,alg_my_sde,alg_spde,alg_jde,
	alg_fep,alg_ifep,
	alg_dgea,
	alg_ede,alg_deeda,alg_dmde,
	alg_eda,
	alg_nsde,
	alg_mosade,
	alg_mode
};// user-specified algorithm type

enum fun_type
{
	sphere=1,
	griewank,
	rastrigin,
	ackey,
	f5,
	rosenbrock,
	step,quartic_with_noise,
	ws_location,
	f2,
	f8,
	camelback,
	f12,
	f13,
	f3,
	f4,
	f14,
	f15,
	foxhole,
	langerman,
	michaelwicz,
	chebyshev,
	zdt3
};// user-specified objective function type

enum stop_type{stoptype_gen=1,stoptype_eval,stoptype_stag,stoptype_delta};
enum ini_type{initype_rnd=1,initype_rng_file,initype_file,initype_ortho};
enum bnd_val_type{bndvaltype_all=1,bndvaltype_file};
enum ob_type{obtype_reflect=1,obtype_absorb,obtype_damp};
enum trunc_type{crowd_dist=1,crowd_entropy,crowd_harmonic};

const int min_SOP_algo_code=1;
const int max_SOP_algo_code=6;
const int min_de_algo_code=7;
const int max_de_algo_code=12;
const int min_ep_algo_code=13;
const int max_ep_algo_code=14;
const int min_ea_algo_code=15;
const int max_ea_algo_code=15;
const int min_hde_algo_code=16;// hybrid de
const int max_hde_algo_code=18;
const int min_eda_algo_code=19;
const int max_eda_algo_code=19;
const int nsde_algo_code=20;

const int mosade_algo_code=21;
const int mode_algo_code=22;
const int min_mode_algo_code=21;
const int max_mode_algo_code=22;
const int mpso_algo_code=2;
const int max_algo_code=22;

inline bool is_pso(int algo_type) { return ((algo_type>=min_SOP_algo_code) && (algo_type<=max_SOP_algo_code)); }
inline bool is_mpso(int algo_type) { return algo_type==mpso_algo_code; }
inline bool is_de(int algo_type) { return ((algo_type>=min_de_algo_code) && (algo_type<=max_de_algo_code)); }
inline bool is_hybrid_de(int algo_type ) { return ((algo_type>=min_hde_algo_code) && (algo_type<=max_hde_algo_code)); }
inline bool is_ep(int algo_type) { return ((algo_type>=min_ep_algo_code) && (algo_type<=max_ep_algo_code)); }
inline bool is_ea(int algo_type) { return ((algo_type>=min_ea_algo_code) && (algo_type<=max_ea_algo_code)); }
inline bool is_eda(int algo_type) { return ((algo_type>=min_eda_algo_code) && (algo_type<=max_eda_algo_code)); }
inline bool is_nsde(int algo_type) { return (nsde_algo_code==algo_type); }
inline bool is_mode(int algo_type) { return ((algo_type>=min_mode_algo_code) && (algo_type<=max_mode_algo_code));  }

inline bool has_max_gen(int stop_type) { return ((stop_type==stoptype_stag) || (stop_type==stoptype_delta)); }

#endif
