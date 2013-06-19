#!/bin/bash

cd pt
echo -n "Testing Poschl-Teller code ... "
../../bin/pt | grep -v '#' > pt.dat
grep -v '#' pt.ref > /tmp/expected
diff pt.dat /tmp/expected
if [ $? -ne 0 ]; then
    echo "failed!"
    diff -y pt.dat /tmp/expected
else
    echo "passed"
fi


