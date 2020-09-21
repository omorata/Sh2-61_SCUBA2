#  lists.mk
#  O. Morata 2020
#
#  rules to list the rules and dependencies in the Makefile used to
#  reduce and analyse the JCMT SCUBA2 data of Sh2-61
#
#  To be included by the main Makefile
#

.SHELLFLAGS = -ec
list:
# lists all the rules in the Makefile
#
	@$(MAKE) -pRrq -f $(MAKEFILE_LIST) : 2> /dev/null |\
            awk -v RS= -F: '/^# File/,/^# Finished Make data base/  \
                {if ($$1 !~ "^[#.]") {print $$1}}' | sort | \
            egrep -v -e '^[^[:alnum:]]' -e '^$@$$' 
##


list_deps:

	@$(MAKE) -pRrq -f $(MAKEFILE_LIST) : 2>/dev/null |\
            awk -v RS= -F: '/^# File/,/^# Finished Make data base/  \
                {if ($$1 !~ "^[#.]") {print $$1" ==> "$$2}}' | sort | \
            egrep -v -e '^[^[:alnum:]]' -e '^$@$$'
##


list_files:
# lists all the possible files created in the rules in the Makefile
#
	@$(MAKE) -pRrq -f $(MAKEFILE_LIST) : 2>/dev/null | \
            awk -v RS= -F: '/^# File/,/^# Finished Make data base/ \
                {if ($$1 !~ "^[#.]") {print $$1}}' | sort | \
            egrep -e '\/'

.PHONY: list list_deps list_files

