# Live progress panel for a pipeline run. Redraws every N seconds.
# Usage: powershell -NoProfile -File tools\show_progress.ps1 -Project "D:\path\to\repo" [-Interval 10]
param(
    [string]$Project = "D:\01_INVESTIGACION\Articulos en curso\IJHG\reproducible_analysis",
    [int]$Interval = 10,
    [string]$Tag = ""
)

if ($Tag) {
    $log = Join-Path $Project "output\logs\$Tag\run.log"
    $modelsDir = Join-Path $Project "output\models\$Tag"
} else {
    $log = Join-Path $Project "output\logs\run.log"
    $modelsDir = Join-Path $Project "output\models"
}

$expected = @(
    "incidence_Wuhan1", "incidence_Wuhan2", "incidence_D614G", "incidence_Beta",
    "incidence_Delta", "incidence_Omicron1", "incidence_Omicron2",
    "mortality_Wuhan1", "mortality_Wuhan2", "mortality_D614G", "mortality_Beta", "mortality_Delta",
    "cfr_Wuhan1", "cfr_Wuhan2", "cfr_D614G", "cfr_Beta", "cfr_Delta"
)

while ($true) {
    Clear-Host
    $now = Get-Date
    Write-Host "  PIPELINE PROGRESS  " -ForegroundColor Black -BackgroundColor Cyan
    Write-Host ("  {0}   {1}" -f $now.ToString("yyyy-MM-dd HH:mm:ss"), (Split-Path $Project -Leaf)) -ForegroundColor DarkGray
    Write-Host ""

    $done = @{}
    if (Test-Path $modelsDir) {
        Get-ChildItem $modelsDir -Filter *.rds -ErrorAction SilentlyContinue |
            ForEach-Object { $done[$_.BaseName] = $_.LastWriteTime }
    }

    $current = "idle"
    if (Test-Path $log) {
        $last = Get-Content $log -Tail 40 -ErrorAction SilentlyContinue
        $fitting = $last | Select-String -Pattern 'fitting (\S+)' | Select-Object -Last 1
        $written = $last | Select-String -Pattern 'models\\(\S+)\.rds' | Select-Object -Last 1
        if ($fitting) {
            $current = $fitting.Matches[0].Groups[1].Value
            if ($written -and $written.LineNumber -gt $fitting.LineNumber) { $current = "between models" }
        }
    }

    $groups = @{ "incidence" = 7; "mortality" = 5; "cfr" = 5 }
    foreach ($oc in @("incidence", "mortality", "cfr")) {
        $set = $expected | Where-Object { $_ -like "${oc}_*" }
        $nDone = ($set | Where-Object { $done.ContainsKey($_) }).Count
        $bar = ""
        foreach ($m in $set) {
            if ($done.ContainsKey($m)) { $bar += "#" } else { $bar += "." }
        }
        $color = if ($nDone -eq $set.Count) { "Green" } else { "Yellow" }
        Write-Host ("  {0,-10} {1,-6} {2}" -f $oc, "$nDone/$($set.Count)", $bar) -ForegroundColor $color
    }

    # ---- parallel workers, if this project runs models through run_parallel.R
    $wdir = Join-Path $Project "output\tables\_workers"
    Write-Host ""
    if (Test-Path $wdir) {
        $running = Test-Path (Join-Path $wdir "RUNNING")
        $stateTxt = if ($running) { "running" } else { "queued, waiting to start" }
        Write-Host ("  PARALLEL  {0}" -f $stateTxt) -ForegroundColor Cyan
        Get-ChildItem $wdir -Directory -ErrorAction SilentlyContinue | Sort-Object Name | ForEach-Object {
            $wd = $_.FullName
            $jobsFile = Join-Path $wd "jobs.txt"
            if (Test-Path $jobsFile) {
                $jobs = @((Get-Content $jobsFile -Raw).Trim() -split ",")
                $nDone = 0
                $now = $null
                foreach ($j in $jobs) {
                    $tag = $j -replace ":", "_"
                    if (Test-Path (Join-Path $wd "$tag.done")) { $nDone++ }
                    elseif (-not $now) { $now = $j }
                }
                $state = if ($nDone -eq $jobs.Count) { "done" }
                         elseif ($now -and $running) { "fitting $now" }
                         elseif ($now) { "next: $now" }
                         else { "starting" }
                $color = if ($nDone -eq $jobs.Count) { "Green" }
                         elseif ($running) { "Yellow" } else { "DarkGray" }
                Write-Host ("    {0,-4} {1}/{2}  {3}" -f $_.Name, $nDone, $jobs.Count, $state) -ForegroundColor $color
            }
        }
    } else {
        Write-Host "  PARALLEL  no plan yet" -ForegroundColor DarkGray
    }

    Write-Host ""
    if ($current -ne "idle" -and $current -ne "between models") {
        Write-Host ("  fitting now: {0}" -f $current) -ForegroundColor Cyan
    } else {
        Write-Host "  no model being fitted (post-processing or idle)" -ForegroundColor Cyan
    }

    $missing = $expected | Where-Object { -not $done.ContainsKey($_) -and $_ -notmatch "Omicron" }
    if ($missing.Count -gt 0 -and $current -eq "between models") {
        Write-Host ("  pending: {0}" -f ($missing -join ", ")) -ForegroundColor DarkYellow
    }

    Write-Host ""
    Write-Host "  last lines" -ForegroundColor DarkGray
    if (Test-Path $log) {
        Get-Content $log -Tail 8 -ErrorAction SilentlyContinue | ForEach-Object {
            Write-Host ("  " + $_) -ForegroundColor DarkGray
        }
    }

    Write-Host ""
    Write-Host ("  refresh every {0}s - Ctrl+C to close" -f $Interval) -ForegroundColor DarkGray
    Start-Sleep -Seconds $Interval
}
