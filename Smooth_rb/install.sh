#!/bin/bash
clear
args=`getopt -o vmpcfib:j: --long make,proper,configure,force,install,build_type:,compiler-cxx:,compiler-c:,out: -n "HDRimage InstallScript" -- "$@"`

if [ $? -ne 0 ];
then
    echo "Error... Usage: $0 -v -m -p -c -f -b[Generation mod]"
    exit 1
fi
eval set -- $args
args=("$@")
CMOD=Release
ROOT_DIRECTORY="$(pwd)"
THREADS=4
OUT_DIR="$(pwd)/bin"

force=false
verbose=""
index=1
for i in $*
do
    case "$i" in
    	--out)
            if [ -n "$i" ]; then
                OUT_DIR=${args[index]}
            fi
          ;;
        -f)
            force=true
          ;;
        -v)
            verbose=1
          ;;
        \?)
            echo -e "\e[00;31mInvalid option\e[00m: -$OPTARG" >&2
            exit 1
          ;;
    esac
    index=$(($index+1))
done
index=1
for i in $*
do
    case "$i" in
        -b)
            if [ -n "$i" ]; then
                CMOD=${args[index]}
            fi
            break
          ;;
    esac
    index=$(($index+1))
done
index=1
CC=gcc
CXX=g++
for i in $*
do
    case "$i" in
        --compiler-c)
            if [ -n "$i" ]; then
                CC=${args[index]}
            fi
            break
          ;;
    esac
    index=$(($index+1))
done
index=1
for i in $*
do
    case "$i" in
        --compiler-cxx)
            if [ -n "$i" ]; then
                CXX=${args[index]}
            fi
            break
          ;;
    esac
    index=$(($index+1))
done
for i in $*
do
    case "$i" in
        -c)
            if $force; then
                echo -e "\e[00;31mRemoving build folder\e[00m"
                cd "$ROOT_DIRECTORY"
                if [ -d "build" ]; then
                    rm -rf build
                fi
            fi
            echo -e "\e[00;31mConfiguring HDRimage\e[00m"
            if [ ! -d "build" ]; then
                mkdir "$ROOT_DIRECTORY/build"
            fi
            cd "$ROOT_DIRECTORY/build"
            env PKG_CONFIG_PATH=~/Dropbox/Yasl/hwloc-1.6.1/build \
            CXX=$CXX CC=$CC cmake -DCMAKE_BUILD_TYPE=$CMOD ..
            break
          ;;
    esac
done
index=1
for i in $*
do
    case "$i" in
        -j)
            if [ -n "$i" ]; then
                THREADS=${args[index]}
            fi
            break
          ;;
    esac
    index=$(($index+1))
done
for i in $*
do
    case "$i" in
        -m|-i)
            echo -e "\e[00;31mGenerating target in $CMOD\e[00m"
            cd "$ROOT_DIRECTORY/build"
            if $force; then
                make clean
            fi
            make -j $THREADS VERBOSE=$verbose
            break
          ;;
    esac
done
for i in $*
do
    case "$i" in
        -i)
            echo -e "\e[00;31mInstalling HDRimage in $OUT_DIR\e[00m"
            cd "$ROOT_DIRECTORY/build"
            make -j $THREADS DESTDIR="$OUT_DIR" install
            break
          ;;
    esac
done
for i in $*
do
    case "$i" in
        -p)
            echo -e "\e[00;31mRemoving build folder\e[00m"
            cd "$ROOT_DIRECTORY"
            if [ -d "build" ]; then
                rm -rf build
            fi
            break
          ;;
    esac
done


cd "$ROOT_DIRECTORY"