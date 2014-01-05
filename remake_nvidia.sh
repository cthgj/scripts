#!/bin/bash

sudo service gdm3 stop 
sudo module-assistant clean nvidia-kernel-source 
sudo module-assistant auto-install nvidia
