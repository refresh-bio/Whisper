all: whisper whisper-index



ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif

ifeq ($(uname_S),Linux)   
	ABI_FLAGS = -fabi-version=6
	ABI_FLAGS_AVX = -fabi-version=0
endif
ifeq ($(uname_S),Darwin)
	ABI_FLAGS = 
	ABI_FLAGS_AVX =
endif


MAPPER_ROOT_DIR = .
MAPPER_MAIN_DIR = mapper
MAPPER_COMMON_DIR = common
MAPPER_PP_DIR = index
MAPPER_LIBS_DIR = libs

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++17 -msse4.2 -mpopcnt $(ABI_FLAGS) -pthread -I $(MAPPER_LIBS_DIR)
CLINK	= -lm -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++17
CFLAGS_AVX = $(CFLAGS) -mavx $(ABI_FLAGS_AVX) -mpopcnt -funroll-loops
CFLAGS_AVX2 = $(CFLAGS) -mavx2 -mfma $(ABI_FLAGS_AVX) -mpopcnt -funroll-loops

# AVX compilation
$(MAPPER_MAIN_DIR)/LevMyers128_AVX.o :  $(MAPPER_MAIN_DIR)/LevMyers128_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/LevMyers128_AVX.cpp -o $@

$(MAPPER_MAIN_DIR)/LevMyers256_AVX.o :  $(MAPPER_MAIN_DIR)/LevMyers256_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/LevMyers256_AVX.cpp -o $@

$(MAPPER_MAIN_DIR)/simd_utils_128_AVX.o :  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX.cpp -o $@	

$(MAPPER_MAIN_DIR)/simd_utils_256_AVX.o :  $(MAPPER_MAIN_DIR)/simd_utils_256_AVX.cpp
	$(CC) $(CFLAGS_AVX) -c  $(MAPPER_MAIN_DIR)/simd_utils_256_AVX.cpp -o $@
	
# AVX2 compilation
$(MAPPER_MAIN_DIR)/LevMyers128_AVX2.o :  $(MAPPER_MAIN_DIR)/LevMyers128_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/LevMyers128_AVX2.cpp -o $@

$(MAPPER_MAIN_DIR)/LevMyers256_AVX2.o :  $(MAPPER_MAIN_DIR)/LevMyers256_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/LevMyers256_AVX2.cpp -o $@

$(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.o :  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.cpp -o $@

$(MAPPER_MAIN_DIR)/simd_utils_256_AVX2.o :  $(MAPPER_MAIN_DIR)/simd_utils_256_AVX2.cpp
	$(CC) $(CFLAGS_AVX2) -c  $(MAPPER_MAIN_DIR)/simd_utils_256_AVX2.cpp -o $@

# general compilation
.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

OBJS := $(MAPPER_MAIN_DIR)/whisper.o \
	$(MAPPER_MAIN_DIR)/bins.o \
	$(MAPPER_MAIN_DIR)/edit_dist.o \
	$(MAPPER_MAIN_DIR)/fastq_reader.o \
	$(MAPPER_MAIN_DIR)/idstore.o \
	$(MAPPER_MAIN_DIR)/indels.o \
	$(MAPPER_MAIN_DIR)/joiner_mgr.o \
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
	$(MAPPER_MAIN_DIR)/params.o \
	$(MAPPER_MAIN_DIR)/queue.o \
	$(MAPPER_MAIN_DIR)/reads.o \
	$(MAPPER_MAIN_DIR)/ref_desc.o \
	$(MAPPER_COMMON_DIR)/reference.o \
	$(MAPPER_MAIN_DIR)/results.o \
	$(MAPPER_MAIN_DIR)/sa.o \
	$(MAPPER_MAIN_DIR)/sam.o \
	$(MAPPER_MAIN_DIR)/sam_part.o \
	$(MAPPER_MAIN_DIR)/sam_matching.o \
	$(MAPPER_MAIN_DIR)/sam_variant.o \
	$(MAPPER_MAIN_DIR)/sam_writer.o \
	$(MAPPER_MAIN_DIR)/simd_utils.o \
	$(MAPPER_MAIN_DIR)/simd_utils_64.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_SSE.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_AVX.o \
	$(MAPPER_MAIN_DIR)/simd_utils_128_AVX2.o \
	$(MAPPER_MAIN_DIR)/simd_utils_256_SSE.o \
	$(MAPPER_MAIN_DIR)/simd_utils_256_AVX.o \
	$(MAPPER_MAIN_DIR)/simd_utils_256_AVX2.o \
	$(MAPPER_MAIN_DIR)/soft_clipping.o \
	$(MAPPER_MAIN_DIR)/stats.o \
	$(MAPPER_COMMON_DIR)/timer.o \
	$(MAPPER_COMMON_DIR)/utils.o \
	$(MAPPER_COMMON_DIR)/variant.o \
	$(MAPPER_MAIN_DIR)/variant_match.o \
	$(MAPPER_LIBS_DIR)/instrset_detect.o
	
INDEX_OBJS := $(MAPPER_PP_DIR)/whisper-index.o \
	$(MAPPER_PP_DIR)/vcf.o \
	$(MAPPER_PP_DIR)/ref_seq.o \
	$(MAPPER_PP_DIR)/suf_arr.o \
	$(MAPPER_PP_DIR)/suf_arr_core.o \
	$(MAPPER_COMMON_DIR)/timer.o \
	$(MAPPER_COMMON_DIR)/utils.o \
	$(MAPPER_COMMON_DIR)/variant.o

whisper: $(OBJS)
	$(CC) $(CLINK) -o $(MAPPER_ROOT_DIR)/$@  $(OBJS) \
	$(MAPPER_LIBS_DIR)/libz.a \
	$(MAPPER_LIBS_DIR)/libbz2.a \
	$(MAPPER_LIBS_DIR)/libaelf64.a \
	$(MAPPER_LIBS_DIR)/libdeflate.a \


whisper-index: $(INDEX_OBJS)
	$(CC) $(CLINK) -o $(MAPPER_ROOT_DIR)/$@  $(INDEX_OBJS)

clean:
	-rm $(MAPPER_MAIN_DIR)/*.o
	-rm $(MAPPER_COMMON_DIR)/*.o
	-rm $(MAPPER_PP_DIR)/*.o
	-rm whisper
	-rm whisper-index
	
