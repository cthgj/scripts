#!/bin/bash
#script to restart networking after failed suspend
echo "sudo modprobe -r e1000e"
sudo modprobe -r e1000e
sleep 3
echo "sudo modprobe e1000e"
sudo modprobe e1000e

