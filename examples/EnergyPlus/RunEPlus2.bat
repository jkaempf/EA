@echo off
: CALL to RUNEPLUS
IF NOT EXIST "%1". MKDIR "%1"
CD %1
IF EXIST in.idf         DEL in.idf
COPY ..\%1.idf .
RunEPlus %1 %2
