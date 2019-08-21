all: whisper whisper-index

MAPPER_ROOT_DIR = .
MAPPER_MAIN_DIR = mapper
MAPPER_COMMON_DIR = common
MAPPER_PP_DIR = index
MAPPER_LIBS_DIR = libs

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++11 -msse2 -fabi-version=6 -pthread -I $(MAPPER_LIBS_DIR)
CLINK	= -lm -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11
CFLAGS_AVX = $(CFLAGS) -mavx -fabi-version=0 -mpopcnt -funroll-loops
CFLAGS_AVX2 = $(CFLAGS) -mavx2 -mfma -fabi-version=0 -mpopcnt -funroll-loops

$(MAPPER_MAIN_DIR)/LevMyers128_AVX.o :  $(MAPPER_MAIN_DIR)/LevMyers128_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/LevMyers128_AVX.cpp -o $@
$(MAPPER_MAIN_DIR)/LevMyers128_AVX2.o :  $(MAPPER_MAIN_DIR)/LevMyers128_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/LevMyers128_AVX2.cpp -o $@
$(MAPPER_MAIN_DIR)/LevMyers256_AVX.o :  $(MAPPER_MAIN_DIR)/LevMyers256_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/LevMyers256_AVX.cpp -o $@
$(MAPPER_MAIN_DIR)/LevMyers256_AVX2.o :  $(MAPPER_MAIN_DIR)/LevMyers256_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/LevMyers256_AVX2.cpp -o $@

$(MAPPER_MAIN_DIR)/simd_utils_128_AVX.o :  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX.cpp -o $@
$(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.o :  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.cpp -o $@
$(MAPPER_MAIN_DIR)/simd_utils_256.o :  $(MAPPER_MAIN_DIR)/simd_utils_256.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/simd_utils_256.cpp -o $@


.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

whisper: $(MAPPER_MAIN_DIR)/whisper.o \
	$(MAPPER_MAIN_DIR)/bins.o \
	$(MAPPER_MAIN_DIR)/edit_dist.o \
	$(MAPPER_COMMON_DIR)/fastq_reader.o \
	$(MAPPER_COMMON_DIR)/idstore.o \
	$(MAPPER_COMMON_DIR)/joiner_mgr.o \
	$(MAPPER_MAIN_DIR)/LevMyers64.o \
	$(MAPPER_MAIN_DIR)/LevMyers128_SSE.o \
	$(MAPPER_MAIN_DIR)/LevMyers128_AVX.o \
	$(MAPPER_MAIN_DIR)/LevMyers128_AVX2.o \
	$(MAPPER_MAIN_DIR)/LevMyers256_AVX.o \
	$(MAPPER_MAIN_DIR)/LevMyers256_AVX2.o \
	$(MAPPER_MAIN_DIR)/mapper.o \
	$(MAPPER_MAIN_DIR)/mapping_core.o \
	$(MAPPER_MAIN_DIR)/mapping_core_sa.o \
	$(MAPPER_MAIN_DIR)/mapping_core_diff_len_reads.o \
	$(MAPPER_MAIN_DIR)/mapping_core_dp.o \
	$(MAPPER_MAIN_DIR)/mapping_core_substrings.o \
	$(MAPPER_MAIN_DIR)/mapping_core_str_comp.o \
	$(MAPPER_COMMON_DIR)/params.o \
	$(MAPPER_COMMON_DIR)/queue.o \
	$(MAPPER_MAIN_DIR)/reads.o \
	$(MAPPER_MAIN_DIR)/ref_desc.o \
	$(MAPPER_COMMON_DIR)/reference.o \
	$(MAPPER_COMMON_DIR)/results.o \
	$(MAPPER_COMMON_DIR)/sa.o \
	$(MAPPER_MAIN_DIR)/sam.o \
	$(MAPPER_MAIN_DIR)/sam_part.o \
	$(MAPPER_MAIN_DIR)/sam_matching.o \
	$(MAPPER_MAIN_DIR)/sam_writer.o \
	$(MAPPER_MAIN_DIR)/simd_utils.o \
	$(MAPPER_MAIN_DIR)/simd_utils_64.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_SSE.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_AVX.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.o \
	$(MAPPER_MAIN_DIR)/simd_utils_256.o \
	$(MAPPER_MAIN_DIR)/soft_clipping.o \
	$(MAPPER_COMMON_DIR)/stats.o \
	$(MAPPER_COMMON_DIR)/timer.o \
	$(MAPPER_COMMON_DIR)/utils.o \
	$(MAPPER_LIBS_DIR)/instrset_detect.o
	$(CC) $(CLINK) -o $(MAPPER_ROOT_DIR)/$@  $(MAPPER_MAIN_DIR)/whisper.o \
	$(MAPPER_MAIN_DIR)/bins.o \
	$(MAPPER_MAIN_DIR)/edit_dist.o \
	$(MAPPER_COMMON_DIR)/fastq_reader.o \
	$(MAPPER_COMMON_DIR)/idstore.o \
	$(MAPPER_COMMON_DIR)/joiner_mgr.o \
	$(MAPPER_MAIN_DIR)/LevMyers64.o \
	$(MAPPER_MAIN_DIR)/LevMyers128_SSE.o \
	$(MAPPER_MAIN_DIR)/LevMyers128_AVX.o \
	$(MAPPER_MAIN_DIR)/LevMyers128_AVX2.o \
	$(MAPPER_MAIN_DIR)/LevMyers256_AVX.o \
	$(MAPPER_MAIN_DIR)/LevMyers256_AVX2.o \
	$(MAPPER_MAIN_DIR)/mapper.o \
	$(MAPPER_MAIN_DIR)/mapping_core.o \
	$(MAPPER_MAIN_DIR)/mapping_core_sa.o \
	$(MAPPER_MAIN_DIR)/mapping_core_diff_len_reads.o \
	$(MAPPER_MAIN_DIR)/mapping_core_dp.o \
	$(MAPPER_MAIN_DIR)/mapping_core_substrings.o \
	$(MAPPER_MAIN_DIR)/mapping_core_str_comp.o \
	$(MAPPER_COMMON_DIR)/params.o \
	$(MAPPER_COMMON_DIR)/queue.o \
	$(MAPPER_MAIN_DIR)/reads.o \
	$(MAPPER_MAIN_DIR)/ref_desc.o \
	$(MAPPER_COMMON_DIR)/reference.o \
	$(MAPPER_COMMON_DIR)/results.o \
	$(MAPPER_COMMON_DIR)/sa.o \
	$(MAPPER_MAIN_DIR)/sam.o \
	$(MAPPER_MAIN_DIR)/sam_part.o \
	$(MAPPER_MAIN_DIR)/sam_matching.o \
	$(MAPPER_MAIN_DIR)/sam_writer.o \
	$(MAPPER_MAIN_DIR)/simd_utils.o \
	$(MAPPER_MAIN_DIR)/simd_utils_64.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_SSE.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_AVX.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.o \
	$(MAPPER_MAIN_DIR)/simd_utils_256.o \
	$(MAPPER_MAIN_DIR)/soft_clipping.o \
	$(MAPPER_COMMON_DIR)/stats.o \
	$(MAPPER_COMMON_DIR)/timer.o \
	$(MAPPER_COMMON_DIR)/utils.o \
	$(MAPPER_LIBS_DIR)/libz.a \
	$(MAPPER_LIBS_DIR)/libbz2.a \
	$(MAPPER_LIBS_DIR)/libaelf64.a \
	$(MAPPER_LIBS_DIR)/libdeflate.a \
	$(MAPPER_LIBS_DIR)/instrset_detect.o

whisper-index: $(MAPPER_PP_DIR)/whisper-index.o \
	$(MAPPER_PP_DIR)/ref_seq.o \
	$(MAPPER_PP_DIR)/suf_arr.o \
	$(MAPPER_PP_DIR)/suf_arr_core.o \
	$(MAPPER_COMMON_DIR)/timer.o \
	$(MAPPER_COMMON_DIR)/utils.o
	$(CC) $(CLINK) -o $(MAPPER_ROOT_DIR)/$@  $(MAPPER_PP_DIR)/whisper-index.o \
	$(MAPPER_PP_DIR)/ref_seq.o \
	$(MAPPER_PP_DIR)/suf_arr.o \
	$(MAPPER_PP_DIR)/suf_arr_core.o \
	$(MAPPER_COMMON_DIR)/timer.o \
	$(MAPPER_COMMON_DIR)/utils.o

clean:
	-rm $(MAPPER_MAIN_DIR)/*.o
	-rm $(MAPPER_COMMON_DIR)/*.o
	-rm $(MAPPER_PP_DIR)/*.o
	-rm whisper
	-rm whisper-index
	
