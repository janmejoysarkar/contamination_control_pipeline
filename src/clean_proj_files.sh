#!/bin/bash
PROJDIR=$(realpath ..)

RAWFD=$PROJDIR/data/raw/full_disk
RAWROI=$PROJDIR/data/raw/roi
PRDFD=$PROJDIR/products/full_disk
PRDROI=$PROJDIR/products/roi

echo $RAWFD
rm -I $RAWFD/*.fits
echo $RAWROI
rm -I $RAWROI/*.fits
echo $PRDFD
rm -I $PRDFD/*.fits
echo $PRDROI
rm -I $PRDROI/*.fits
