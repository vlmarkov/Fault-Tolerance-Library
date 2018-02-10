#!/bin/bash

./makeit.sh && cd build && qsub torque.job && qsub torque-kill.job
