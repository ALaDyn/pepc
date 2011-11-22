<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html"/>
  <xsl:template match="/benchrun">
    <html><head>
    <link rel="stylesheet" type="text/css" href="style.css"/>
    <title>Report</title></head>
    <body bgcolor="#EEEEEE">
      <h1><a name="topOfPage">JuBE Report</a></h1>
      <xsl:apply-templates select="jobstartdate"/>
      <xsl:apply-templates select="jobenddate"/>

      <!-- table witch some general information -->
      <br></br>
      <h2>General Information</h2>
      <table>
        <tr>
          <td><a href="#identifier"><xsl:text>Identifier</xsl:text></a></td>
          <td><xsl:value-of select="benchmark/@identifier"/></td>
        </tr>
        <tr>
          <td><a href="#ncpus"><xsl:text>Number Of CPUs</xsl:text></a></td>
          <td><xsl:value-of select="benchmark/@ncpus"/></td>
        </tr>
        <tr>
          <td><a href="#tasks"><xsl:text>Number Of Tasks</xsl:text></a></td>
          <td><xsl:value-of select="compile/params/@tasks"/></td>
        </tr>
        <tr>
          <td><a href="#taskspernode"><xsl:text>Tasks Per Node</xsl:text></a></td>
          <td><xsl:value-of select="compile/params/@taskspernode"/></td>
        </tr>
        <tr>
          <td><a href="#threadspertask"><xsl:text>Threads Per Task</xsl:text></a></td>
          <td><xsl:value-of select="compile/params/@threadspertask"/></td>
        </tr>
        <tr>
          <td><a href="#walltime"><xsl:text>Wall Clock Time</xsl:text></a></td>
          <td><xsl:value-of select="analyse/walltime/@value"/><xsl:text> sec.</xsl:text></td>
        </tr>

      </table>



      <xsl:apply-templates select="benchmark"/> 
      <xsl:apply-templates select="compile/params"/>
      <xsl:apply-templates select="analyse"/>
      <xsl:apply-templates select="stdoutfile"/>
      <xsl:apply-templates select="stderrfile"/>
      <!--  <xsl:apply-templates/> -->
    </body ></html>
  </xsl:template>
  
  <xsl:template match="compile/params">
    <h2>Compile Section</h2>
    <p><a href="#topOfPage">Top Of Page</a></p>
    <xsl:text> Compile command:</xsl:text>
    <xsl:value-of select="../command"/>
    <br></br>
    <div style="width:100%;">
      <table width="100%" border="2" cellspacing="0" cellpadding="0" style="overflow:hidden; table-layout:fixed;word-wrap:break-word;">
        <tr style="font-weight:bold">
          <th>Parameter</th>
          <th>Value</th>
        </tr>
        <xsl:for-each select="@*">
          <tr>
            <xsl:choose>
              <xsl:when test="boolean(name()='tasks')">
                <td><a name="tasks"><xsl:value-of select="name()" /></a></td>
              </xsl:when>
              <xsl:when test="boolean(name()='taskspernode')">
                <td><a name="taskspernode"><xsl:value-of select="name()" /></a></td>
              </xsl:when>
              <xsl:when test="boolean(name()='threadspertask')">
                <td><a name="threadspertask"><xsl:value-of select="name()" /></a></td>
              </xsl:when>
              <xsl:otherwise>       
              <td><xsl:value-of select="name()" /></td>
            </xsl:otherwise>
          </xsl:choose>

          <xsl:choose>
            <xsl:when test="boolean(name()='execname')">
              <!-- shorten string for execname -->
              <xsl:variable name="execname_short1">
                <xsl:value-of select="substring-after(current(),'tmp/')"/>
              </xsl:variable>
              <xsl:variable name="execname_short2">
                <xsl:value-of select="substring-after($execname_short1,'/')"/>
              </xsl:variable>
              <xsl:variable name="execname_veryshort">
                <xsl:value-of select="substring-after($execname_short2,'/')"/>
              </xsl:variable>
              
              <td><xsl:value-of select="$execname_veryshort"/></td>
            </xsl:when>
            <xsl:when test="boolean(name()='outdir')">
              <td><xsl:value-of select="substring-after(current(),'tmp/')"/></td>
            </xsl:when>
            <xsl:when test="boolean(name()='subdir')">
              <td><xsl:value-of select="substring-after(current(),'tmp/')"/></td>
            </xsl:when>
            <xsl:when test="boolean(name()='rundir')">
              <td><xsl:value-of select="substring-after(current(),'tmp/')"/></td>
            </xsl:when>
            <xsl:otherwise>
              <td><xsl:value-of select="." /></td>
            </xsl:otherwise>
          </xsl:choose>
        </tr>       
      </xsl:for-each>
    </table>
  </div>
  <br></br>
</xsl:template>

<xsl:template match="bstartdate">
  <h2>XXX started: </h2>
  <xsl:value-of select="."/>
  <br></br>
</xsl:template>

<xsl:template match="benddate">
  <h2>XXX ended: </h2>
  <xsl:value-of select="."/>
  <br></br>
</xsl:template> 
  
<xsl:template match="analyse">
  <br></br>
  <h2>Analyse Section</h2>
  <p><a href="#topOfPage">Top Of Page</a></p>
  <table border="2" width="100%">
    <tr style="font-weight:bold">
      <th>Name</th>
      <th>Value</th>
      <th>Unit</th>
    </tr>    
    
    <xsl:for-each select="*">
      <tr>
        <xsl:choose>
          <xsl:when test="boolean(name()='walltime')">
            <td><a name="walltime"><xsl:value-of select="name()"/></a></td>
          </xsl:when>
          <xsl:otherwise>
            <td><xsl:value-of select="name()"/></td>
          </xsl:otherwise>
        </xsl:choose>
        <xsl:choose>
          <xsl:when test="boolean(values)">
            <td>
              <table width="100%">
                <tr>
                  <th>Method</th>
                  <th>Value</th>
                </tr>
                <xsl:for-each select="values">
                  <tr>
                    <td class="multiValueElement"><xsl:value-of select="@name"/></td>
                    <td class="multiValueElement"><xsl:value-of select="@value"/></td>
                  </tr>
                </xsl:for-each>
              </table>
            </td>
            <td></td>
            </xsl:when>
            <xsl:otherwise>
              
              
              <td><xsl:value-of select="@value"/></td>
              <td><xsl:value-of select="@unit"/></td>
            </xsl:otherwise>
          </xsl:choose>
        </tr>
      </xsl:for-each>
  </table>
  <br></br>



  <!-- end second part -->
  <br></br>
</xsl:template>
  
<xsl:template match="benchmark">
  <br></br>
  <h2>Benchrun</h2>
  <p><a href="#topOfPage">Top Of Page</a></p>
  <table  width="100%" border="2" cellspacing="0" cellpadding="0" style="overflow:hidden; table-layout:fixed;">
    <tr style="font-weight:bold">
      <th>Name</th>
      <th>Value</th>
    </tr>
    <xsl:for-each select="@*">
      <tr>
        <xsl:choose>
          <xsl:when test="boolean(name()='ncpus')">
            <td> <a name="ncpus"><xsl:value-of select="name()" /></a></td>
          </xsl:when>
          <xsl:when test="boolean(name()='identifier')">
            <td> <a name="identifier"><xsl:value-of select="name()" /></a></td>
          </xsl:when>
          <xsl:otherwise>
            <td> <xsl:value-of select="name()" /></td>
          </xsl:otherwise>
        </xsl:choose>
        
        <xsl:choose>          
          <xsl:when test="boolean(name()='stderrfile')">
            <td><xsl:value-of select="substring-after(current(),'logs/')"/></td>
          </xsl:when>
          <xsl:when test="boolean(name()='stdoutfile')">
            <td><xsl:value-of select="substring-after(current(),'logs/')"/></td>
          </xsl:when>
          <xsl:when test="boolean(name()='subdir')">
            <td><xsl:value-of select="substring-after(current(),'tmp/')"/></td>
          </xsl:when>
          <xsl:otherwise>
            <td><xsl:value-of select="." /></td>
          </xsl:otherwise>
        </xsl:choose>
        
      </tr>
    </xsl:for-each>
  </table>
  <br></br>
</xsl:template>

<xsl:template match="stdoutfile">
  <h2>Standard Output File:</h2>
  <p><a href="#topOfPage">Top Of Page</a></p>
  <textarea name="output" cols="150" rows="30">
    <xsl:value-of select="."/>
  </textarea>
</xsl:template>

<xsl:template match="stderrfile">
  <h2>Standard Error File</h2>
  <p><a href="#topOfPage">Top Of Page</a></p>
  <textarea name="output" cols="150" rows="30">
    <xsl:value-of select="."/>
  </textarea>
</xsl:template>

<xsl:template match="jobstartdate">
  <h2>Calculation started: <xsl:value-of select="."/></h2>
</xsl:template>

<xsl:template match="jobenddate">
  <h2>Calculation ended: <xsl:value-of select="."/></h2>
</xsl:template>
</xsl:stylesheet>
