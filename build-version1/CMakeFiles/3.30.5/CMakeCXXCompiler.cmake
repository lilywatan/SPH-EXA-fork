set(CMAKE_CXX_COMPILER "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/bin/g++")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "GNU")
set(CMAKE_CXX_COMPILER_VERSION "13.2.0")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "17")
set(CMAKE_CXX_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_CXX_STANDARD_LATEST "23")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters;cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates;cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates;cxx_std_17;cxx_std_20;cxx_std_23")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates")
set(CMAKE_CXX17_COMPILE_FEATURES "cxx_std_17")
set(CMAKE_CXX20_COMPILE_FEATURES "cxx_std_20")
set(CMAKE_CXX23_COMPILE_FEATURES "cxx_std_23")
set(CMAKE_CXX26_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_CXX_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_CXX_COMPILER_AR "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/bin/gcc-ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_LINKER_LINK "")
set(CMAKE_LINKER_LLD "")
set(CMAKE_CXX_COMPILER_LINKER "/usr/bin/ld")
set(CMAKE_CXX_COMPILER_LINKER_ID "GNU")
set(CMAKE_CXX_COMPILER_LINKER_VERSION 2.34)
set(CMAKE_CXX_COMPILER_LINKER_FRONTEND_VARIANT GNU)
set(CMAKE_MT "")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_COMPILER_IS_GNUCXX 1)
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;mpp;CPP;ixx;cppm;ccm;cxxm;c++m)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)

foreach (lang IN ITEMS C OBJC OBJCXX)
  if (CMAKE_${lang}_COMPILER_ID_RUN)
    foreach(extension IN LISTS CMAKE_${lang}_SOURCE_FILE_EXTENSIONS)
      list(REMOVE_ITEM CMAKE_CXX_SOURCE_FILE_EXTENSIONS ${extension})
    endforeach()
  endif()
endforeach()

set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)
set(CMAKE_CXX_LINKER_DEPFILE_SUPPORTED FALSE)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/ncurses-6.5-6ps2tp4jwxsy7mlj2euw7tnbqqudwsv3/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/curl-8.10.1-xp7lugzcaezdnq7r6gtzk3l2mormtufb/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openssl-3.4.0-krtm3zj5zzkvmdfi4ep4s6gt7cnrddgu/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/zlib-ng-2.2.1-doufymgdkq7brialss4zbxorxljbc2ib/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/nghttp2-1.63.0-6i2325xcljmb36ocvkojn53rstism33n/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/hdf5-1.12.2-qz66wcdd6apl5snokh4y6isvjobvp4eq/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/pkgconf-1.9.5-4e2bo3gvfuivabceric2dxmk4fk4cciy/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-4.1.3-efaqdqjshm5pznecwekvbmt42qbwhinl/include;/cluster/23-11-2-1/slurm-23-11-2-1/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/pmix-4.1.2-zndjodtzfji4gyda5nnm3o3cec4c3pfd/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libedit-3.1-20210216-ispc5jpqugoxdmjqve44idnrwmx3zxh7/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/numactl-2.0.14-hfslryfsz5m4e7ukdn2a4oq3wkrhzcot/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libevent-2.1.12-hvi4ccst4bkgqh5k4hhqfpnepp2enajm/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/hwloc-2.7.1-eeyokp3qomnvy73qrfhaduyxdyq2orzu/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libxml2-2.9.13-ndomtw6tqe4plooe4xpjd43e3sowopsr/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/zlib-1.2.12-j4b6zegwseutr44qyr66ym767esbxvjm/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/xz-5.2.5-mhrz5subhwqlf35mzkrrigebxkkb7bje/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libiconv-1.16-pdflaobqhkm2yizzmiscm3g2hpqnlowy/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libpciaccess-0.16-6i5pp5phelvy7kdi633adayyabn62lkp/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/python-3.13.0-liiwqoibau2rrs6p3tvznmfxbp5xsuu5/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/util-linux-uuid-2.40.2-qnvnyd5rqo5gczrdiwrlyuklpqkmclvu/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/sqlite-3.46.0-gyjmfhriwkazs7r4zx5nhowbjttuw67w/include;/cluster/24-05-1-1/libffi-3.4.6/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gettext-0.22.5-szntn2mmu3qsupcosibh7p37w5udpoq5/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/zstd-1.5.6-73jijrqavjggtlrxht4d3iframk575g3/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gdbm-1.23-ypb4vc6yush2sphgowqy23uyx5h764hi/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/readline-8.2-y2kxjpiddslljrcxixx6dkluqsdlees7/include;/cluster/24-05-1-1/expat-2.6.2/include;/cluster/24-05-1-1/libbsd-0.12.2/include;/cluster/24-05-1-1/libmd-1.0.4/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/bzip2-1.0.8-x6snebgh4jj5rwwlh7dbhq7n5fapgevj/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mpc-1.3.1-wffm27cmf5mug7uaweowp3s47twjgxub/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mpfr-4.2.0-qoadsiirmvnb5thwkub72c4zxvebf6rg/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gmp-6.2.1-fggj7dy7hu7kstyq62znpd67tizyyjza/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/include/c++/13.2.0;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/include/c++/13.2.0/x86_64-pc-linux-gnu;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/include/c++/13.2.0/backward;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/lib/gcc/x86_64-pc-linux-gnu/13.2.0/include;/usr/local/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/include;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/lib/gcc/x86_64-pc-linux-gnu/13.2.0/include-fixed/x86_64-linux-gnu;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/lib/gcc/x86_64-pc-linux-gnu/13.2.0/include-fixed;/usr/include/x86_64-linux-gnu;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openssl-3.4.0-krtm3zj5zzkvmdfi4ep4s6gt7cnrddgu/lib64;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/lib64;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/lib/gcc/x86_64-pc-linux-gnu/13.2.0;/lib/x86_64-linux-gnu;/lib64;/usr/lib/x86_64-linux-gnu;/usr/lib64;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/ncurses-6.5-6ps2tp4jwxsy7mlj2euw7tnbqqudwsv3/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/curl-8.10.1-xp7lugzcaezdnq7r6gtzk3l2mormtufb/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/zlib-ng-2.2.1-doufymgdkq7brialss4zbxorxljbc2ib/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/nghttp2-1.63.0-6i2325xcljmb36ocvkojn53rstism33n/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/hdf5-1.12.2-qz66wcdd6apl5snokh4y6isvjobvp4eq/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/pkgconf-1.9.5-4e2bo3gvfuivabceric2dxmk4fk4cciy/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-4.1.3-efaqdqjshm5pznecwekvbmt42qbwhinl/lib;/cluster/23-11-2-1/slurm-23-11-2-1/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/pmix-4.1.2-zndjodtzfji4gyda5nnm3o3cec4c3pfd/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libedit-3.1-20210216-ispc5jpqugoxdmjqve44idnrwmx3zxh7/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/numactl-2.0.14-hfslryfsz5m4e7ukdn2a4oq3wkrhzcot/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libevent-2.1.12-hvi4ccst4bkgqh5k4hhqfpnepp2enajm/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/hwloc-2.7.1-eeyokp3qomnvy73qrfhaduyxdyq2orzu/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libxml2-2.9.13-ndomtw6tqe4plooe4xpjd43e3sowopsr/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/zlib-1.2.12-j4b6zegwseutr44qyr66ym767esbxvjm/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/xz-5.2.5-mhrz5subhwqlf35mzkrrigebxkkb7bje/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libiconv-1.16-pdflaobqhkm2yizzmiscm3g2hpqnlowy/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/libpciaccess-0.16-6i5pp5phelvy7kdi633adayyabn62lkp/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/python-3.13.0-liiwqoibau2rrs6p3tvznmfxbp5xsuu5/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/util-linux-uuid-2.40.2-qnvnyd5rqo5gczrdiwrlyuklpqkmclvu/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/sqlite-3.46.0-gyjmfhriwkazs7r4zx5nhowbjttuw67w/lib;/cluster/24-05-1-1/libffi-3.4.6/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gettext-0.22.5-szntn2mmu3qsupcosibh7p37w5udpoq5/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/zstd-1.5.6-73jijrqavjggtlrxht4d3iframk575g3/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gdbm-1.23-ypb4vc6yush2sphgowqy23uyx5h764hi/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/readline-8.2-y2kxjpiddslljrcxixx6dkluqsdlees7/lib;/cluster/24-05-1-1/expat-2.6.2/lib;/cluster/24-05-1-1/libbsd-0.12.2/lib;/cluster/24-05-1-1/libmd-1.0.4/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/bzip2-1.0.8-x6snebgh4jj5rwwlh7dbhq7n5fapgevj/lib;/cluster/24-05-1-1/gcc-runtime-9.3.0/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gcc-13.2.0-hgptpx2eoraipaxlrxijwyj5jxznibqq/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mpc-1.3.1-wffm27cmf5mug7uaweowp3s47twjgxub/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mpfr-4.2.0-qoadsiirmvnb5thwkub72c4zxvebf6rg/lib;/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/gmp-6.2.1-fggj7dy7hu7kstyq62znpd67tizyyjza/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
set(CMAKE_CXX_COMPILER_CLANG_RESOURCE_DIR "")

set(CMAKE_CXX_COMPILER_IMPORT_STD "")
### Imported target for C++23 standard library
set(CMAKE_CXX23_COMPILER_IMPORT_STD_NOT_FOUND_MESSAGE "Unsupported generator: Unix Makefiles")



