#!/usr/bin/env bash

xterm -e '/bin/bash --rcfile ~/1.rc' && wmctrl -r :ACTIVE: -b add,maximized_vert &
xterm -e '/bin/bash --rcfile ~/2.rc' &
xterm -e '/bin/bash --rcfile ~/3.rc' &
xterm -e '/bin/bash --rcfile ~/4.rc' &
