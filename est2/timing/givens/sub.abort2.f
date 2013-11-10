      subroutine abort2

* Purpose:  To take the place of calls to ABORT, which is not present
*           on some computers.  This makes the code more portable.

      stop 'Abnormal termination. See standard output for error message'
      end
