#!/usr/bin/env bash
/home/eric/Projects/sampler/tmp/cmake-build-debug/MTmgen cg_model_1.txt
echo "... model done..."
/home/eric/Projects/sampler/tmp/cmake-build-debug/MTrgen -C cg_model_1.cfg -c on
echo "... simulator configuration done..."
/home/eric/Projects/sampler/tmp/cmake-build-debug/MTrgen -C cg_model_1.cfg
echo "... simulation done..."
evince cg_model_1_rep.pdf
echo "... all done."