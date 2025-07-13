# stdbuf -oL bash run-convzx-conv-tests.sh |& tee foo.txt

cat $0

exe=~/tmp/winstall/bin/woodland_squirrel_unittest
ls -ltrh $exe

nthr=1

tnres=4; rnres=6

function run {
    echo $cmp $opts
    time OMP_NUM_THREADS=$nthr $exe $cmp $opts
}

cmp=ct_we_zx
for hs in 0 1; do
    for i in 20 21 22; do
        opts="testcase=$i,halfspace=$hs,srfrecon=0,nres=$tnres"; run
    done
done

hs=0
cmp=ct_we_zx
for i in 20 21 22; do    
    opts="testcase=$i,halfspace=$hs,flatelem=1,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,flatelem=1,dislocorder=0,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,flatelem=1,ntri=4,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,flatelem=1,ntri=4,dislocorder=0,nres=$tnres"; run

    opts="testcase=$i,halfspace=$hs,dislocorder=0,srfrecon=0,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,dislocorder=1,srfrecon=0,nres=$tnres"; run

    opts="testcase=$i,halfspace=$hs,srfrecon=0,ntri=4,nres=$tnres"; run

    opts="testcase=$i,halfspace=$hs,exacttan=1,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,tanorder=2,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,dislocorder=1,tanorder=2,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,tanorder=4,nres=$tnres"; run
    opts="testcase=$i,halfspace=$hs,c2spline=1,nres=$tnres"; run
done

cmp=ct_oe_zx
for hs in 0 1; do
    for i in 20 21 22 10 11 12; do
        opts="testcase=$i,halfspace=$hs,nres=$rnres"; run
    done
done

cmp=ct_we_zx
for hs in 0 1; do
    for i in 10 11 12; do
        opts="testcase=$i,halfspace=$hs,srfrecon=0,nres=$tnres"; run
        opts="testcase=$i,halfspace=$hs,srfrecon=0,dislocorder=0,nres=$tnres"; run
    done
done
