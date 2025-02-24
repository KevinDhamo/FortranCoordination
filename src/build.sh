#!/bin/bash

# Function to print usage
print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Build the coordination analysis program"
    echo ""
    echo "Options:"
    echo "  -h, --help        Show this help message"
    echo "  -j, --jobs N      Use N parallel jobs for compilation (default: number of CPU cores)"
    echo "  --debug           Build with debug flags"
    echo "  --profile         Build with profiling enabled"
    echo "  --cluster         Build with maximum optimization for cluster"
    echo "  --clean          Clean build artifacts before building"
    echo ""
    echo "Example:"
    echo "  $0 --cluster -j 8    # Build cluster version using 8 cores"
}

# Function to detect CPU capabilities
detect_cpu() {
    if command -v lscpu >/dev/null 2>&1; then
        CPU_MODEL=$(lscpu | grep "Model name" | cut -d: -f2 | xargs)
        CPU_CORES=$(lscpu | grep "^CPU(s):" | cut -d: -f2 | xargs)
        if lscpu | grep -q "avx512"; then
            HAS_AVX512=1
            echo "AVX-512 support detected"
        else
            HAS_AVX512=0
            echo "No AVX-512 support found"
        fi
    else
        CPU_MODEL="Unknown"
        CPU_CORES=$(nproc 2>/dev/null || echo 1)
        HAS_AVX512=0
    fi
}

# Function to detect system capabilities
detect_system() {
    # Detect memory
    if [ -f /proc/meminfo ]; then
        TOTAL_RAM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        TOTAL_RAM_GB=$((TOTAL_RAM_KB / 1024 / 1024))
    else
        TOTAL_RAM_GB=1
    fi

    # Detect if running on a cluster
    if [ $CPU_CORES -gt 16 ] || [ $TOTAL_RAM_GB -gt 64 ]; then
        IS_CLUSTER=1
        echo "System appears to be a cluster node"
    else
        IS_CLUSTER=0
        echo "System appears to be a standard workstation"
    fi
}

# Function to check requirements
check_requirements() {
    local missing_req=0

    # Check for gfortran
    if ! command -v gfortran >/dev/null 2>&1; then
        echo "Error: gfortran is not installed"
        echo "Install with: sudo apt-get install gfortran"
        missing_req=1
    fi

    # Check for make
    if ! command -v make >/dev/null 2>&1; then
        echo "Error: make is not installed"
        echo "Install with: sudo apt-get install make"
        missing_req=1
    fi

    # Check OpenMP support
    if ! gfortran -v 2>&1 | grep -q "enable-libgomp"; then
        echo "Warning: OpenMP support might not be available"
    fi

    if [ $missing_req -eq 1 ]; then
        exit 1
    fi
}

# Parse command line arguments
BUILD_TYPE="release"
PARALLEL_JOBS=0
CLEAN_FIRST=0

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            print_usage
            exit 0
            ;;
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --debug)
            BUILD_TYPE="debug"
            shift
            ;;
        --profile)
            BUILD_TYPE="profile"
            shift
            ;;
        --cluster)
            BUILD_TYPE="cluster"
            shift
            ;;
        --clean)
            CLEAN_FIRST=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Main script
echo "=== Coordination Analysis Build Script ==="

# Check requirements
check_requirements

# Detect CPU and system capabilities
detect_cpu
detect_system

echo "System Detection Results:"
echo "  CPU Model: $CPU_MODEL"
echo "  CPU Cores: $CPU_CORES"
echo "  Total RAM: $TOTAL_RAM_GB GB"
echo "  AVX-512 Support: $HAS_AVX512"

# Set number of parallel jobs if not specified
if [ $PARALLEL_JOBS -eq 0 ]; then
    PARALLEL_JOBS=$CPU_CORES
fi

# Clean if requested
if [ $CLEAN_FIRST -eq 1 ]; then
    echo "Cleaning previous build..."
    make clean
fi

# Determine optimal build type if not explicitly specified
if [ "$BUILD_TYPE" = "release" ] && [ $IS_CLUSTER -eq 1 ]; then
    echo "Automatically selecting cluster build type based on system detection"
    BUILD_TYPE="cluster"
fi

# Build the program
echo "Building with type: $BUILD_TYPE"
echo "Using $PARALLEL_JOBS parallel jobs"

# Select appropriate make target based on CPU capabilities
case $BUILD_TYPE in
    cluster)
        if [ $HAS_AVX512 -eq 1 ]; then
            make cluster -j$PARALLEL_JOBS
            EXECUTABLE="coordination_analysis_cluster"
        else
            make cluster_old -j$PARALLEL_JOBS
            EXECUTABLE="coordination_analysis_cluster"
        fi
        ;;
    debug)
        make debug -j$PARALLEL_JOBS
        EXECUTABLE="coordination_analysis_debug"
        ;;
    profile)
        make profile -j$PARALLEL_JOBS
        EXECUTABLE="coordination_analysis_prof"
        ;;
    *)
        make release -j$PARALLEL_JOBS
        EXECUTABLE="coordination_analysis"
        ;;
esac

# Check if build was successful
if [ $? -eq 0 ] && [ -f "$EXECUTABLE" ]; then
    echo "Build successful! Executable: $EXECUTABLE"
    echo ""
    echo "To run the program:"
    echo "  ./$EXECUTABLE"
    
    # Print additional information for cluster builds
    if [ "$BUILD_TYPE" = "cluster" ]; then
        echo ""
        echo "Cluster optimization information:"
        echo "  - Built with advanced vectorization support"
        if [ $HAS_AVX512 -eq 1 ]; then
            echo "  - Using AVX-512 instructions"
        fi
        echo "  - Optimized for parallel execution"
        echo "  - Recommended MPI ranks: $CPU_CORES"
    fi
else
    echo "Build failed!"
    exit 1
fi