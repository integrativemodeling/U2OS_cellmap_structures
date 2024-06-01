from IMP.bbm.restraints import parameterize

parameter.params_from_inferred_binding_modes(
    pdb_files, out_prefix,
    chains=("A", "B"),
    dbscan_eps=0.1,
    dbscan_min_samples=5,
    lr=0.0001, lrq=0.0001, 
    check_tol=0.001, conv_tol=0.000001,
    n_m_steps=1, n_steps=100,
    **kwargs
)
