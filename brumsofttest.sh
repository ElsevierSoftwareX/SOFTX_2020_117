#!/bin/sh

while getopts s:u:t:g: option
do
 case "${option}"
 in
 s) SERVER=${OPTARG};;
 u) USER=${OPTARG};;
 t) TEST=${OPTARG};;
 g) GITDATA=$OPTARG;;
 esac
done

SUBURL="${SERVER}/cli/submit_test/user/${USER}/testProfile/${TEST}/gitcheckout/${GITDATA}"
echo "Submit URL: "$SUBURL

SESSIONID="`wget -qO- $SUBURL`"
echo "SessionID: "$SESSIONID

STATUSURL="${SERVER}/cli/get_status_only/session/${SESSIONID}/user/${USER}"
DATAURL="${SERVER}/cli/get_status/session/${SESSIONID}/user/${USER}"

echo "Status URL: "$STATUSURL
echo "Data URL: "$DATAURL

EXIT=0

while [ $EXIT -eq 0 ]; do
  STATUS="`wget -qO- $STATUSURL`"
  if [ "$STATUS" = "Passed" ]; then
    echo $STATUS
    echo "======"
    DATA="`wget -qO- $DATAURL`"
    echo $DATA
    exit 0
  fi
  if [ "$STATUS" = "Failed" ]; then
    echo $STATUS
    echo "======"
    DATA="`wget -qO- $DATAURL`"
    echo $DATA
    exit 1
  fi
  sleep 0.1
done