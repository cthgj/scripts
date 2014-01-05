#!/bin/bash

killall -9 plasma-desktop;
kbuildsycoca4 --noincremental;
plasma-desktop;
