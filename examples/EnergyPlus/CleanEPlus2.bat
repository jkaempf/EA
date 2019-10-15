@echo off
: CLEAN after RUNEPLUS
COPY .\%1\%1.eso .
RMDIR /s /q %1
