# Debug Data Flow Issue

## Problem
The app is processing 15,357 genes instead of the uploaded file data (which should only have 5 genes from test_mouse_data.csv).

## Likely Causes
1. **Browser Session Cache**: Previous app session data is cached
2. **Shiny Reactive Values**: Old data persists in reactive values
3. **Test Data Override**: Automatic test data is overriding uploaded data
4. **File Upload Not Working**: The upload process isn't replacing the data properly

## Debugging Steps

### Step 1: Clear Browser Cache
1. Close all browser tabs with the app
2. Clear browser cache completely
3. Restart browser
4. Navigate to app again

### Step 2: Use Clear Data Button
1. In development mode, click "Clear All Data" button
2. Verify data is cleared
3. Upload test_mouse_data.csv
4. Check console output

### Step 3: Check Console Debug Output
After uploading test_mouse_data.csv, you should see:
```
🔍 DEBUG: Processed matrix stored with dimensions: 5 x 6
🔍 DEBUG: Sample genes in processed data: ENSMUSG00000020108, ENSMUSG00000027490, ENSMUSG00000025746
```

Then in DESeq2 analysis:
```
🔍 DEBUG: Expression data dimensions: 5 x 6
🔍 DEBUG: Sample gene IDs: ENSMUSG00000020108, ENSMUSG00000027490, ENSMUSG00000025746
```

### Step 4: If Still Seeing 15,357 Genes
The issue is that the reactive values are not being properly updated. This could be:
1. The file upload is failing silently
2. Some other process is overriding the data
3. The app is using cached data

## Expected vs Actual

### Expected (test_mouse_data.csv):
- 5 genes: ENSMUSG00000020108, ENSMUSG00000027490, ENSMUSG00000025746, ENSMUSG00000020125, ENSMUSG00000033845
- 6 samples: Sample1, Sample2, Sample3, Sample4, Sample5, Sample6

### If You See 15,357 Genes:
This is definitely NOT from your uploaded file. The test mouse data only has 5 genes.

## Immediate Solution
1. Use "Clear All Data" button in development mode
2. Refresh the page completely (Ctrl+F5 or Cmd+Shift+R)
3. Upload test_mouse_data.csv again
4. Check the debug output in console

## Debug Output to Look For
```
🔍 DEBUG: Processed matrix stored with dimensions: 5 x 6
🔍 DEBUG: Sample genes in processed data: ENSMUSG00000020108, ENSMUSG00000027490, ENSMUSG00000025746
🔍 DEBUG: Expression data dimensions: 5 x 6  
🔍 DEBUG: Sample gene IDs: ENSMUSG00000020108, ENSMUSG00000027490, ENSMUSG00000025746
🔍 Detected mouse gene IDs (5 ENSMUSG patterns)
```

If you don't see this debug output, the file upload is not working properly.