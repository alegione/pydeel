#!/bin/bash

set -e
errors=0

# Run unit tests
python pydeel/pydeel_test.py || {
    echo "'python python/pydeel/pydeel_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E pydeel/*.py || {
    echo 'pylint -E pydeel/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
