/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_xelem_call.h
 *  timestamp: 2009-12-03 09:12:56 +0100
 *  
 */


#define xelem_name                                       keys
#define xelem_name_packed                                key
#define xelem_type_c                                     key_type_c
#define xelem_size_c                                     1
#define xelem_type_mpi                                   key_type_mpi
#define xelem_size_mpi                                   key_size_mpi
#define xelem_sltype_t                                   slkey_t
#define xelem_byte                                       key_byte
#define xelem_assign(_s_, _d_)                           key_assign((_s_)->keys, (_d_)->keys)
#define xelem_assign_at(_s_, _sat_, _d_)                 key_assign_at((_s_)->keys, _sat_, (_d_)->keys)
#define xelem_null(_e_)                                  key_null((_e_)->keys)
#define xelem_inc(_e_)                                   key_inc((_e_)->keys)
#define xelem_dec(_e_)                                   key_dec((_e_)->keys)
#define xelem_add(_e_, _n_)                              key_add((_e_)->keys, _n_)
#define xelem_sub(_e_, _n_)                              key_sub((_e_)->keys, _n_)
#define xelem_copy(_s_, _d_)                             key_copy((_s_)->keys, (_d_)->keys)
#define xelem_ncopy(_s_, _d_, _n_)                       key_ncopy((_s_)->keys, (_d_)->keys, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       key_nmove((_s_)->keys, (_d_)->keys, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            key_copy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      key_ncopy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      key_nmove_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   key_xchange((_e0_)->keys, (_e1_)->keys, (_t_)->keys)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  key_xchange_at((_e0_)->keys, _at0_, (_e1_)->keys, _at1_, (_t_)->keys)
#define xelem_buf(_e_)                                   (_e_)->keys
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->keys[_at_])
#define xelem_mpi_datatype                               pepckeys_key_mpi_datatype

#if !defined(xelem_key_not)
# ifdef xelem_call_key
xelem_call_key
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_key_not
xelem_call_key_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       indices
#define xelem_name_packed                                index
#define xelem_type_c                                     index_type_c
#define xelem_size_c                                     1
#define xelem_type_mpi                                   index_type_mpi
#define xelem_size_mpi                                   index_size_mpi
#define xelem_sltype_t                                   slindex_t
#define xelem_byte                                       index_byte
#define xelem_assign(_s_, _d_)                           index_assign((_s_)->indices, (_d_)->indices)
#define xelem_assign_at(_s_, _sat_, _d_)                 index_assign_at((_s_)->indices, _sat_, (_d_)->indices)
#define xelem_null(_e_)                                  index_null((_e_)->indices)
#define xelem_inc(_e_)                                   index_inc((_e_)->indices)
#define xelem_dec(_e_)                                   index_dec((_e_)->indices)
#define xelem_add(_e_, _n_)                              index_add((_e_)->indices, _n_)
#define xelem_sub(_e_, _n_)                              index_sub((_e_)->indices, _n_)
#define xelem_copy(_s_, _d_)                             index_copy((_s_)->indices, (_d_)->indices)
#define xelem_ncopy(_s_, _d_, _n_)                       index_ncopy((_s_)->indices, (_d_)->indices, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       index_nmove((_s_)->indices, (_d_)->indices, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            index_copy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      index_nmove_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   index_xchange((_e0_)->indices, (_e1_)->indices, (_t_)->indices)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  index_xchange_at((_e0_)->indices, _at0_, (_e1_)->indices, _at1_, (_t_)->indices)
#define xelem_buf(_e_)                                   (_e_)->indices
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->indices[_at_])
#define xelem_mpi_datatype                               pepckeys_index_mpi_datatype

#if defined(SL_INDEX) && !defined(xelem_index_not)
# ifdef xelem_call_index
xelem_call_index
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_index_not
xelem_call_index_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


/* DATAX_TEMPLATE_BEGIN */
#define xelem_name                                       data0
#define xelem_name_packed                                data0
#define xelem_type_c                                     data0_type_c
#define xelem_size_c                                     data0_size_c
#define xelem_type_mpi                                   data0_type_mpi
#define xelem_size_mpi                                   data0_size_mpi
#define xelem_sltype_t                                   sldata0_t
#define xelem_byte                                       data0_byte
#define xelem_assign(_s_, _d_)                           data0_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data0_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data0_null(_e_)
#define xelem_inc(_e_)                                   data0_inc(_e_)
#define xelem_dec(_e_)                                   data0_dec(_e_)
#define xelem_add(_e_, _n_)                              data0_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data0_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data0_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data0_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data0_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data0_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data0_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data0_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data0_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data0_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data0
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data0[(_at_) * data0_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[0]

#if defined(SL_DATA0) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data1
#define xelem_name_packed                                data1
#define xelem_type_c                                     data1_type_c
#define xelem_size_c                                     data1_size_c
#define xelem_type_mpi                                   data1_type_mpi
#define xelem_size_mpi                                   data1_size_mpi
#define xelem_sltype_t                                   sldata1_t
#define xelem_byte                                       data1_byte
#define xelem_assign(_s_, _d_)                           data1_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data1_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data1_null(_e_)
#define xelem_inc(_e_)                                   data1_inc(_e_)
#define xelem_dec(_e_)                                   data1_dec(_e_)
#define xelem_add(_e_, _n_)                              data1_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data1_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data1_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data1_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data1_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data1_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data1_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data1_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data1_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data1_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data1
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data1[(_at_) * data1_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[1]

#if defined(SL_DATA1) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data2
#define xelem_name_packed                                data2
#define xelem_type_c                                     data2_type_c
#define xelem_size_c                                     data2_size_c
#define xelem_type_mpi                                   data2_type_mpi
#define xelem_size_mpi                                   data2_size_mpi
#define xelem_sltype_t                                   sldata2_t
#define xelem_byte                                       data2_byte
#define xelem_assign(_s_, _d_)                           data2_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data2_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data2_null(_e_)
#define xelem_inc(_e_)                                   data2_inc(_e_)
#define xelem_dec(_e_)                                   data2_dec(_e_)
#define xelem_add(_e_, _n_)                              data2_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data2_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data2_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data2_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data2_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data2_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data2_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data2_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data2_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data2_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data2
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data2[(_at_) * data2_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[2]

#if defined(SL_DATA2) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data3
#define xelem_name_packed                                data3
#define xelem_type_c                                     data3_type_c
#define xelem_size_c                                     data3_size_c
#define xelem_type_mpi                                   data3_type_mpi
#define xelem_size_mpi                                   data3_size_mpi
#define xelem_sltype_t                                   sldata3_t
#define xelem_byte                                       data3_byte
#define xelem_assign(_s_, _d_)                           data3_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data3_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data3_null(_e_)
#define xelem_inc(_e_)                                   data3_inc(_e_)
#define xelem_dec(_e_)                                   data3_dec(_e_)
#define xelem_add(_e_, _n_)                              data3_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data3_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data3_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data3_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data3_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data3_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data3_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data3_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data3_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data3_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data3
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data3[(_at_) * data3_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[3]

#if defined(SL_DATA3) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data4
#define xelem_name_packed                                data4
#define xelem_type_c                                     data4_type_c
#define xelem_size_c                                     data4_size_c
#define xelem_type_mpi                                   data4_type_mpi
#define xelem_size_mpi                                   data4_size_mpi
#define xelem_sltype_t                                   sldata4_t
#define xelem_byte                                       data4_byte
#define xelem_assign(_s_, _d_)                           data4_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data4_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data4_null(_e_)
#define xelem_inc(_e_)                                   data4_inc(_e_)
#define xelem_dec(_e_)                                   data4_dec(_e_)
#define xelem_add(_e_, _n_)                              data4_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data4_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data4_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data4_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data4_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data4_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data4_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data4_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data4_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data4_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data4
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data4[(_at_) * data4_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[4]

#if defined(SL_DATA4) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data5
#define xelem_name_packed                                data5
#define xelem_type_c                                     data5_type_c
#define xelem_size_c                                     data5_size_c
#define xelem_type_mpi                                   data5_type_mpi
#define xelem_size_mpi                                   data5_size_mpi
#define xelem_sltype_t                                   sldata5_t
#define xelem_byte                                       data5_byte
#define xelem_assign(_s_, _d_)                           data5_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data5_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data5_null(_e_)
#define xelem_inc(_e_)                                   data5_inc(_e_)
#define xelem_dec(_e_)                                   data5_dec(_e_)
#define xelem_add(_e_, _n_)                              data5_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data5_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data5_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data5_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data5_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data5_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data5_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data5_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data5_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data5_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data5
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data5[(_at_) * data5_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[5]

#if defined(SL_DATA5) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data6
#define xelem_name_packed                                data6
#define xelem_type_c                                     data6_type_c
#define xelem_size_c                                     data6_size_c
#define xelem_type_mpi                                   data6_type_mpi
#define xelem_size_mpi                                   data6_size_mpi
#define xelem_sltype_t                                   sldata6_t
#define xelem_byte                                       data6_byte
#define xelem_assign(_s_, _d_)                           data6_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data6_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data6_null(_e_)
#define xelem_inc(_e_)                                   data6_inc(_e_)
#define xelem_dec(_e_)                                   data6_dec(_e_)
#define xelem_add(_e_, _n_)                              data6_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data6_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data6_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data6_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data6_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data6_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data6_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data6_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data6_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data6_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data6
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data6[(_at_) * data6_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[6]

#if defined(SL_DATA6) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data7
#define xelem_name_packed                                data7
#define xelem_type_c                                     data7_type_c
#define xelem_size_c                                     data7_size_c
#define xelem_type_mpi                                   data7_type_mpi
#define xelem_size_mpi                                   data7_size_mpi
#define xelem_sltype_t                                   sldata7_t
#define xelem_byte                                       data7_byte
#define xelem_assign(_s_, _d_)                           data7_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data7_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data7_null(_e_)
#define xelem_inc(_e_)                                   data7_inc(_e_)
#define xelem_dec(_e_)                                   data7_dec(_e_)
#define xelem_add(_e_, _n_)                              data7_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data7_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data7_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data7_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data7_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data7_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data7_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data7_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data7_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data7_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data7
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data7[(_at_) * data7_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[7]

#if defined(SL_DATA7) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data8
#define xelem_name_packed                                data8
#define xelem_type_c                                     data8_type_c
#define xelem_size_c                                     data8_size_c
#define xelem_type_mpi                                   data8_type_mpi
#define xelem_size_mpi                                   data8_size_mpi
#define xelem_sltype_t                                   sldata8_t
#define xelem_byte                                       data8_byte
#define xelem_assign(_s_, _d_)                           data8_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data8_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data8_null(_e_)
#define xelem_inc(_e_)                                   data8_inc(_e_)
#define xelem_dec(_e_)                                   data8_dec(_e_)
#define xelem_add(_e_, _n_)                              data8_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data8_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data8_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data8_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data8_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data8_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data8_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data8_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data8_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data8_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data8
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data8[(_at_) * data8_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[8]

#if defined(SL_DATA8) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data9
#define xelem_name_packed                                data9
#define xelem_type_c                                     data9_type_c
#define xelem_size_c                                     data9_size_c
#define xelem_type_mpi                                   data9_type_mpi
#define xelem_size_mpi                                   data9_size_mpi
#define xelem_sltype_t                                   sldata9_t
#define xelem_byte                                       data9_byte
#define xelem_assign(_s_, _d_)                           data9_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data9_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data9_null(_e_)
#define xelem_inc(_e_)                                   data9_inc(_e_)
#define xelem_dec(_e_)                                   data9_dec(_e_)
#define xelem_add(_e_, _n_)                              data9_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data9_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data9_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data9_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data9_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data9_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data9_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data9_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data9_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data9_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data9
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data9[(_at_) * data9_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[9]

#if defined(SL_DATA9) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data10
#define xelem_name_packed                                data10
#define xelem_type_c                                     data10_type_c
#define xelem_size_c                                     data10_size_c
#define xelem_type_mpi                                   data10_type_mpi
#define xelem_size_mpi                                   data10_size_mpi
#define xelem_sltype_t                                   sldata10_t
#define xelem_byte                                       data10_byte
#define xelem_assign(_s_, _d_)                           data10_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data10_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data10_null(_e_)
#define xelem_inc(_e_)                                   data10_inc(_e_)
#define xelem_dec(_e_)                                   data10_dec(_e_)
#define xelem_add(_e_, _n_)                              data10_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data10_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data10_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data10_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data10_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data10_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data10_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data10_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data10_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data10_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data10
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data10[(_at_) * data10_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[10]

#if defined(SL_DATA10) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data11
#define xelem_name_packed                                data11
#define xelem_type_c                                     data11_type_c
#define xelem_size_c                                     data11_size_c
#define xelem_type_mpi                                   data11_type_mpi
#define xelem_size_mpi                                   data11_size_mpi
#define xelem_sltype_t                                   sldata11_t
#define xelem_byte                                       data11_byte
#define xelem_assign(_s_, _d_)                           data11_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data11_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data11_null(_e_)
#define xelem_inc(_e_)                                   data11_inc(_e_)
#define xelem_dec(_e_)                                   data11_dec(_e_)
#define xelem_add(_e_, _n_)                              data11_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data11_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data11_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data11_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data11_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data11_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data11_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data11_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data11_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data11_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data11
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data11[(_at_) * data11_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[11]

#if defined(SL_DATA11) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data12
#define xelem_name_packed                                data12
#define xelem_type_c                                     data12_type_c
#define xelem_size_c                                     data12_size_c
#define xelem_type_mpi                                   data12_type_mpi
#define xelem_size_mpi                                   data12_size_mpi
#define xelem_sltype_t                                   sldata12_t
#define xelem_byte                                       data12_byte
#define xelem_assign(_s_, _d_)                           data12_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data12_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data12_null(_e_)
#define xelem_inc(_e_)                                   data12_inc(_e_)
#define xelem_dec(_e_)                                   data12_dec(_e_)
#define xelem_add(_e_, _n_)                              data12_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data12_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data12_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data12_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data12_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data12_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data12_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data12_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data12_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data12_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data12
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data12[(_at_) * data12_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[12]

#if defined(SL_DATA12) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data13
#define xelem_name_packed                                data13
#define xelem_type_c                                     data13_type_c
#define xelem_size_c                                     data13_size_c
#define xelem_type_mpi                                   data13_type_mpi
#define xelem_size_mpi                                   data13_size_mpi
#define xelem_sltype_t                                   sldata13_t
#define xelem_byte                                       data13_byte
#define xelem_assign(_s_, _d_)                           data13_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data13_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data13_null(_e_)
#define xelem_inc(_e_)                                   data13_inc(_e_)
#define xelem_dec(_e_)                                   data13_dec(_e_)
#define xelem_add(_e_, _n_)                              data13_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data13_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data13_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data13_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data13_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data13_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data13_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data13_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data13_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data13_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data13
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data13[(_at_) * data13_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[13]

#if defined(SL_DATA13) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data14
#define xelem_name_packed                                data14
#define xelem_type_c                                     data14_type_c
#define xelem_size_c                                     data14_size_c
#define xelem_type_mpi                                   data14_type_mpi
#define xelem_size_mpi                                   data14_size_mpi
#define xelem_sltype_t                                   sldata14_t
#define xelem_byte                                       data14_byte
#define xelem_assign(_s_, _d_)                           data14_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data14_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data14_null(_e_)
#define xelem_inc(_e_)                                   data14_inc(_e_)
#define xelem_dec(_e_)                                   data14_dec(_e_)
#define xelem_add(_e_, _n_)                              data14_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data14_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data14_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data14_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data14_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data14_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data14_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data14_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data14_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data14_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data14
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data14[(_at_) * data14_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[14]

#if defined(SL_DATA14) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data15
#define xelem_name_packed                                data15
#define xelem_type_c                                     data15_type_c
#define xelem_size_c                                     data15_size_c
#define xelem_type_mpi                                   data15_type_mpi
#define xelem_size_mpi                                   data15_size_mpi
#define xelem_sltype_t                                   sldata15_t
#define xelem_byte                                       data15_byte
#define xelem_assign(_s_, _d_)                           data15_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data15_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data15_null(_e_)
#define xelem_inc(_e_)                                   data15_inc(_e_)
#define xelem_dec(_e_)                                   data15_dec(_e_)
#define xelem_add(_e_, _n_)                              data15_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data15_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data15_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data15_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data15_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data15_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data15_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data15_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data15_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data15_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data15
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data15[(_at_) * data15_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[15]

#if defined(SL_DATA15) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data16
#define xelem_name_packed                                data16
#define xelem_type_c                                     data16_type_c
#define xelem_size_c                                     data16_size_c
#define xelem_type_mpi                                   data16_type_mpi
#define xelem_size_mpi                                   data16_size_mpi
#define xelem_sltype_t                                   sldata16_t
#define xelem_byte                                       data16_byte
#define xelem_assign(_s_, _d_)                           data16_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data16_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data16_null(_e_)
#define xelem_inc(_e_)                                   data16_inc(_e_)
#define xelem_dec(_e_)                                   data16_dec(_e_)
#define xelem_add(_e_, _n_)                              data16_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data16_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data16_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data16_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data16_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data16_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data16_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data16_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data16_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data16_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data16
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data16[(_at_) * data16_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[16]

#if defined(SL_DATA16) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data17
#define xelem_name_packed                                data17
#define xelem_type_c                                     data17_type_c
#define xelem_size_c                                     data17_size_c
#define xelem_type_mpi                                   data17_type_mpi
#define xelem_size_mpi                                   data17_size_mpi
#define xelem_sltype_t                                   sldata17_t
#define xelem_byte                                       data17_byte
#define xelem_assign(_s_, _d_)                           data17_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data17_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data17_null(_e_)
#define xelem_inc(_e_)                                   data17_inc(_e_)
#define xelem_dec(_e_)                                   data17_dec(_e_)
#define xelem_add(_e_, _n_)                              data17_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data17_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data17_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data17_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data17_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data17_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data17_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data17_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data17_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data17_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data17
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data17[(_at_) * data17_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[17]

#if defined(SL_DATA17) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data18
#define xelem_name_packed                                data18
#define xelem_type_c                                     data18_type_c
#define xelem_size_c                                     data18_size_c
#define xelem_type_mpi                                   data18_type_mpi
#define xelem_size_mpi                                   data18_size_mpi
#define xelem_sltype_t                                   sldata18_t
#define xelem_byte                                       data18_byte
#define xelem_assign(_s_, _d_)                           data18_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data18_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data18_null(_e_)
#define xelem_inc(_e_)                                   data18_inc(_e_)
#define xelem_dec(_e_)                                   data18_dec(_e_)
#define xelem_add(_e_, _n_)                              data18_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data18_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data18_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data18_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data18_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data18_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data18_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data18_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data18_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data18_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data18
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data18[(_at_) * data18_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[18]

#if defined(SL_DATA18) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
#define xelem_name                                       data19
#define xelem_name_packed                                data19
#define xelem_type_c                                     data19_type_c
#define xelem_size_c                                     data19_size_c
#define xelem_type_mpi                                   data19_type_mpi
#define xelem_size_mpi                                   data19_size_mpi
#define xelem_sltype_t                                   sldata19_t
#define xelem_byte                                       data19_byte
#define xelem_assign(_s_, _d_)                           data19_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data19_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data19_null(_e_)
#define xelem_inc(_e_)                                   data19_inc(_e_)
#define xelem_dec(_e_)                                   data19_dec(_e_)
#define xelem_add(_e_, _n_)                              data19_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data19_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data19_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data19_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data19_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data19_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data19_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data19_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data19_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data19_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_buf(_e_)                                   (_e_)->data19
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data19[(_at_) * data19_size_c])
#define xelem_mpi_datatype                               pepckeys_data_mpi_datatype[19]

#if defined(SL_DATA19) && !defined(xelem_data_not)
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype
/* DATAX_TEMPLATE_END */

#undef xelem_key_not
#undef xelem_index_not
#undef xelem_data_not
#undef xelem_call_key
#undef xelem_call_key_not
#undef xelem_call_index
#undef xelem_call_index_not
#undef xelem_call_data
#undef xelem_call_data_not
#undef xelem_call
#undef xelem_call_not
