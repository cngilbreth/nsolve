#!/bin/bash

cd ho
echo -n "Testing harmonic oscillator code ... "
../../bin/ho | grep -v '#' > ho.dat
grep -v '#' ho.ref > /tmp/expected
diff ho.dat /tmp/expected
if [ $? -ne 0 ]; then
    echo "failed!"
    diff -y ho.dat /tmp/expected
    exit 1
else
    echo "passed"
fi


