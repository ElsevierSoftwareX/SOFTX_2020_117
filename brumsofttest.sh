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

case $SESSIONID in
    ''|*[!0-9]*) exit 1 ;;
    *) echo "ID Okay" ;;
esac

STATUSURL="${SERVER}/cli/get_status_only/session/${SESSIONID}/user/${USER}"
DATAURL="${SERVER}/cli/get_status/session/${SESSIONID}/user/${USER}"

echo "Status URL: "$STATUSURL
echo "Data URL: "$DATAURL

EXIT=0

while [ $EXIT -eq 0 ]; do
  STATUS="`wget -qO- $STATUSURL`"
  case "$STATUS" in
    "Passed")
      echo '\n\n'$STATUS
      echo "======"
      DATA="`wget -qO- $DATAURL`"
      echo $DATA
      exit 0
      ;;
    "Failed")
      echo '\n\n'$STATUS
      echo "======"
      DATA="`wget -qO- $DATAURL`"
      echo $DATA
      exit 1
      ;;
    "Queued") printf "*" ;;
    "Running") printf "." ;;
    *)
      echo "Unknown response: "$DATA
      exit 1
      ;;
  esac
  sleep 1
done
