/* Dark mode toggle — clicking the nav brand icon (#nav-brand-icon) toggles
 * body.dark-mode, persists preference to localStorage, and inverts the
 * background mouse images + the icon itself via CSS. */
$(document).ready(function () {
  var DARK_KEY = 'meDashboard_darkMode';

  // Restore preference on load
  if (localStorage.getItem(DARK_KEY) === '1') {
    $('body').addClass('dark-mode');
  }

  // Click on the mouse logo beside "mouseExperiment" to toggle dark mode
  $(document).on('click', '#nav-brand-icon', function () {
    var isDark = $('body').toggleClass('dark-mode').hasClass('dark-mode');
    localStorage.setItem(DARK_KEY, isDark ? '1' : '0');
  });
});
