#!/bin/sh

sage --preparse qform.sage
mv qform.sage.py qform.py
sage main.sage
